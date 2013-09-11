/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SliceImage.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SliceImage.h"

#include <ctime>
#include <cmath>
#include <limits>

#include "vtkLine.h"
#include "vtkDataArray.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLineSource.h"
#include "vtkObjectFactory.h"
#include "vtkParametricSpline.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtksys/ios/sstream>

vtkStandardNewMacro(SliceImage);

void SliceImage::SetCurve(vtkPolyDataAlgorithm *line_segment) {
	this->Spline = NULL;
	this->LineSegment = line_segment;
	m_clip_type = LINE;

	this->Modified();
}
void SliceImage::SetCurve(vtkParametricSpline *spline) {
	this->Spline = spline;
	//this->LineSegment = NULL;

	m_clip_type = SPLINE;
	this->Modified();
}

//----------------------------------------------------------------------------
// Construct an instance of SliceImage filter.
SliceImage::SliceImage()
{
	this->Spline = NULL;
	this->LineSegment = NULL;

	// by default process active point scalars
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
void SliceImage::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int SliceImage::RequestInformation(vtkInformation* information,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
	using namespace std;
	// Get input and output pipeline information.
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// Get the input whole extent.
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), m_input_extent);

	inInfo->Get(vtkDataObject::ORIGIN(), m_input_origin);
	inInfo->Get(vtkDataObject::SPACING(), m_input_spacing);

	m_n_pixels = 1;
	std::vector<SplinePoint> piecewise_linear_curve = LinearizeClippingFunction();
	if (!piecewise_linear_curve.empty())
		m_n_pixels = piecewise_linear_curve.back().cumm_length / m_input_spacing[0];

	m_output_extent[0] = 0;
	m_output_extent[1] = m_n_pixels;
	m_output_extent[2] = m_input_extent[4];
	m_output_extent[3] = m_input_extent[5];
	m_output_extent[4] = 0;
	m_output_extent[5] = 0;

	// Store the new whole extent for the output.
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), m_output_extent, 6);
	information->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), m_output_extent, 6);

	// Set the number of point data componets to the number of
	// components in the gradient vector.
	vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);

	return 1;
}

//----------------------------------------------------------------------------
// This method computes the input extent necessary to generate the output.
int SliceImage::RequestUpdateExtent(vtkInformation*,
		vtkInformationVector** inputVector,
		vtkInformationVector* /*outputVector*/)
{
	// NGHK: Optimize region
	// Get input and output pipeline information.
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// Get the input whole extent.
	int wholeExtent[6];
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);

	// Store the update extent needed from the intput.
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), wholeExtent, 6);

	return 1;
}

void SliceImage::get(int delta, int &x, int &y) const
{
	float fdelta = delta*1.0/m_n_pixels;
	if (fdelta < 0) fdelta = 0.0f;
	if (fdelta > 1) fdelta = 1.0f;
	if (m_pt0_in_pixels[0] == m_pt1_in_pixels[0])
		x = m_pt0_in_pixels[0];
	else
		x = m_pt0_in_pixels[0] + (m_pt1_in_pixels[0]-m_pt0_in_pixels[0]) * fdelta;

	if (m_pt0_in_pixels[1] == m_pt1_in_pixels[1])
		y = m_pt0_in_pixels[1];
	else
		y = m_pt0_in_pixels[1] + (m_pt1_in_pixels[1]-m_pt0_in_pixels[1]) * fdelta;
}
bool SliceImage::WorldToPixelCoordinates(const double world_pos[3], int pixel_pos[3])
{
	bool in_domain = true;
	for (int i=0; i<3; ++i) {
		pixel_pos[i] = std::floor((world_pos[i]-m_input_origin[i]) / m_input_spacing[i] + 0.5);

		if (pixel_pos[i] < m_input_extent[2*i]) {
			pixel_pos[i] = m_input_extent[2*i];
			in_domain = false;
		} else if (pixel_pos[i] > m_input_extent[2*i+1]) {
			pixel_pos[i] = m_input_extent[2*i+1];
			in_domain = false;
		}
	}

	return in_domain;
}

//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values
// out of extent.
template <class T>
void SliceImage::SliceImageExecute(SliceImage *self,
		vtkImageData *inData, T *inPtr,
		vtkImageData *outData, double *outPtr,
		int outExt[6], int id)
{
	int maxX, maxY, maxZ;
	vtkIdType inIncX, inIncY, inIncZ;
	vtkIdType outIncX, outIncY, outIncZ;
	unsigned long count = 0;
	unsigned long target;
	int *inExt = inData->GetExtent();
	//int *wholeExtent;
	vtkIdType *inIncs;
	double r[3];
	//int useZMin, useZMax, useYMin, useYMax, useXMin, useXMax;

	// find the region to loop over
	maxX = outExt[1] - outExt[0];
	maxY = outExt[3] - outExt[2];
	maxZ = outExt[5] - outExt[4];
	target = static_cast<unsigned long>((maxZ+1)*(maxY+1)/50.0);
	target++;

	// Get increments to march through data
	inData->GetContinuousIncrements(outExt, inIncX, inIncY, inIncZ);
	outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

	// The data spacing is important for computing the gradient.
	// central differences (2 * ratio).
	// Negative because below we have (min - max) for dx ...
	inData->GetSpacing(r);
	r[0] = -0.5 / r[0];
	r[1] = -0.5 / r[1];
	r[2] = -0.5 / r[2];

	// get some other info we need
	inIncs = inData->GetIncrements();
	//wholeExtent = outData->GetExtent();

	// Move the pointer to the correct starting position.
	inPtr += (outExt[0]-inExt[0])*inIncs[0] +
			(outExt[2]-inExt[2])*inIncs[1] +
			(outExt[4]-inExt[4])*inIncs[2];

	// Loop through output pixels
	int x_indices[maxX+1];
	std::vector<SplinePoint> spline_pts = LinearizeClippingFunction();
	double world_pt[3];
	double delta_length = spline_pts.back().cumm_length / maxX;
	size_t spline_idx = 1;
	for (int x=0; x<=maxX; x++) {
		int pixel_coord[3];
		while ((spline_idx+1 < spline_pts.size()) && (spline_pts[spline_idx].cumm_length < x*delta_length)) {
			spline_idx ++;
		}

		// Interpolate linearly
		double alpha =
				(x*delta_length - spline_pts[spline_idx-1].cumm_length) /
				(spline_pts[spline_idx].cumm_length - spline_pts[spline_idx-1].cumm_length);
		for (int dir=0; dir<3; ++dir) {
			world_pt[dir] =
					(1-alpha) * spline_pts[spline_idx-1].pt[dir] + alpha * spline_pts[spline_idx].pt[dir];
		}

		if (this->WorldToPixelCoordinates(world_pt, pixel_coord)) {
			x_indices[x] = inIncs[0]*pixel_coord[0] + inIncs[1]*pixel_coord[1];
		} else
			x_indices[x] = -1;
	}

	for (int z=0; z<=maxZ; z++)
	{
		for (int y=0; y<=maxY; y++)
		{
			if (!id)
			{
				if (!(count%target)) self->UpdateProgress(count/(50.0*target));
				count++;
			}

			int z_offset = y * inIncs[2];

			for (int x=0; x<=maxX; x++)
			{
				if (x_indices[x] >= 0)
					*outPtr = inPtr[z_offset + x_indices[x]];
				else
					*outPtr = std::numeric_limits<double>::quiet_NaN();
				++outPtr;
			}
			outPtr += outIncY;
		}
		outPtr += outIncZ;
	}
}

int SliceImage::RequestData(
		vtkInformation* request,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
	if (!this->Superclass::RequestData(request, inputVector, outputVector))
		return 0;

	vtkImageData* output = vtkImageData::GetData(outputVector);
	vtkDataArray* outArray = output->GetPointData()->GetScalars();
	outArray->SetName((outArray->GetName()?outArray->GetName():"value"));
	return 1;
}

//----------------------------------------------------------------------------
// This method contains a switch statement that calls the correct
// templated function for the input data type.  This method does handle
// boundary conditions.
void SliceImage::ThreadedRequestData(vtkInformation*,
		vtkInformationVector** inputVector,
		vtkInformationVector*,
		vtkImageData*** inData,
		vtkImageData** outData,
		int outExt[6],
		int threadId)
{
	// Get the input and output data objects.
	vtkImageData* input = inData[0][0];
	vtkImageData* output = outData[0];

	// The ouptut scalar type must be double to store proper gradients.
	if(output->GetScalarType() != VTK_DOUBLE)
	{
		vtkErrorMacro("Execute: output ScalarType is "
				<< output->GetScalarType() << "but must be double.");
		return;
	}

	vtkDataArray* inputArray = this->GetInputArrayToProcess(0, inputVector);
	if (!inputArray)
	{
		vtkErrorMacro("No input array was found. Cannot execute");
		return;
	}

	if(inputArray->GetNumberOfComponents() != 1)
	{
		vtkErrorMacro(
				"Execute: input has more than one component.");
		return;
	}

	void* inPtr = inputArray->GetVoidPointer(0);
	double* outPtr = static_cast<double *>(
			output->GetScalarPointerForExtent(outExt));
	switch(inputArray->GetDataType())
	{
	vtkTemplateMacro(
			SliceImageExecute(this, input, static_cast<VTK_TT*>(inPtr),
					output, outPtr, outExt, threadId)
	);
	default:
		vtkErrorMacro("Execute: Unknown ScalarType " << input->GetScalarType());
		return;
	}
}

unsigned long SliceImage::GetMTime()
{
	unsigned long time = vtkThreadedImageAlgorithm::GetMTime();
	if (LineSegment != NULL)
		time = std::max(time, LineSegment->GetMTime());
	if (Spline != NULL)
		time = std::max(time, Spline->GetMTime());
	return time;
}

std::vector<SliceImage::SplinePoint> SliceImage::LinearizeClippingFunction() {
	std::vector<SplinePoint> result;
	switch (m_clip_type) {
	case LINE: {
		vtkLineSource *line = vtkLineSource::SafeDownCast(this->LineSegment);
		if (line == NULL)
			return result;

		SplinePoint pt0, pt1;
		line->GetPoint1(pt0.pt);
		line->GetPoint2(pt1.pt);

		result.push_back(pt0);
		result.push_back(pt1);
		break;
	}
	case SPLINE: {
		if (Spline == NULL)
			return result;

		double u[3], Du[9];
		u[0] = u[1] = u[2] = 0.0;

		// The 32 is an arbitrary sampling density
		size_t n_pts = 2;
		if (Spline->GetPoints()->GetNumberOfPoints() >= 2)
			n_pts = 32*(Spline->GetPoints()->GetNumberOfPoints()-1); // > 2

		double dt = 1.0/(n_pts-1);

		result.resize(n_pts);

		for (size_t i=0; i<n_pts; ++i) {
			u[0] = i * dt;
			Spline->Evaluate(u, result[i].pt, Du);
		}
		break;
	}
	}

	for (size_t i=0; i<result.size(); ++i) result[i].pt[2] = 0.0;

	result[0].cumm_length = 0;

	for (size_t i=1; i<result.size(); ++i) {
		double sq_length = 0.0;
		for (int j=0; j<3; ++j)
			sq_length += (result[i-1].pt[j]-result[i].pt[j])*(result[i-1].pt[j]-result[i].pt[j]);

		result[i].cumm_length = result[i-1].cumm_length + sqrt(sq_length);
	}

	return result;
}
