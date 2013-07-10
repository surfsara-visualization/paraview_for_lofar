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

#include "vtkLine.h"
#include "vtkDataArray.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLineSource.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtksys/ios/sstream>

vtkStandardNewMacro(SliceImage);
vtkCxxSetObjectMacro(SliceImage, CutFunction, vtkImplicitFunction);
vtkCxxSetObjectMacro(SliceImage, Curve, vtkPolyDataAlgorithm);

//----------------------------------------------------------------------------
// Construct an instance of SliceImage filter.
SliceImage::SliceImage()
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	this->CutFunction = NULL;
	this->Curve = NULL;

	// by default process active point scalars
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
void SliceImage::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int SliceImage::RequestInformation(vtkInformation*,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
  using namespace std;
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	// Get input and output pipeline information.
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// Get the input whole extent.
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), m_extent);

	vtkLineSource *line = vtkLineSource::SafeDownCast(this->GetCurve());
	if (line != NULL)
	{
		{
			double input_origin[3];
			inInfo->Get(vtkDataObject::ORIGIN(), input_origin);
			double *spacing = inInfo->Get(vtkDataObject::SPACING());

			double *pt0 = line->GetPoint1();
			double *pt1 = line->GetPoint2();

			for (int i=0; i<3; ++i)
			{
				std::cout << "pt[" << i << "]: " << pt0[i] << " .. " << pt1[i] << std::endl;

				m_pt0_in_pixels[i] = m_extent[2*i] + (pt0[i]-input_origin[i]) / spacing[i];
				m_pt0_in_pixels[i] = std::min(std::max(m_pt0_in_pixels[i], m_extent[2*i]), m_extent[2*i+1]);

				m_pt1_in_pixels[i] = m_extent[2*i] + (pt1[i]-input_origin[i]) / spacing[i];
				m_pt1_in_pixels[i] = std::min(std::max(m_pt1_in_pixels[i], m_extent[2*i]), m_extent[2*i+1]);
			}

			std::cout << "pt0_in_pixels: " << m_pt0_in_pixels[0] << ", " << m_pt0_in_pixels[1] << ", " << m_pt0_in_pixels[2] << std::endl;
			std::cout << "pt1_in_pixels: " << m_pt1_in_pixels[0] << ", " << m_pt1_in_pixels[1] << ", " << m_pt1_in_pixels[2] << std::endl;
		}

		m_n_pixels = std::max(
				abs(m_pt1_in_pixels[0]-m_pt0_in_pixels[0]),
				abs(m_pt1_in_pixels[1]-m_pt0_in_pixels[1]));

		m_extent[0] = 0;
		m_extent[1] = m_n_pixels;
		m_extent[2] = m_extent[4];
		m_extent[3] = m_extent[5];
		m_extent[4] = 0;
		m_extent[5] = 0;

		if (m_n_pixels == 0)
			for (int i=0; i<6; ++i)
				m_extent[i] = 0;

	} else {
		m_extent[0] = 0;
		m_extent[1] = 42;
		m_extent[2] = 0;
		m_extent[3] = 42;
		m_extent[4] = 0;
		m_extent[5] = 42;
	}

	// Store the new whole extent for the output.
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), m_extent, 6);

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

//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values
// out of extent.
template <class T>
void SliceImageExecute(SliceImage *self,
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

	// Loop through ouput pixels
	int x_indices[maxX+1];
	for (int x=0; x<=maxX; x++) {
		int inX, inY;
		self->get(x+outExt[0], inX, inY);
		x_indices[x] = inX * inIncs[0] + inY * inIncs[1];
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
				*outPtr = inPtr[z_offset + x_indices[x]];
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
	{
		return 0;
	}
	vtkImageData* output = vtkImageData::GetData(outputVector);
	vtkDataArray* outArray = output->GetPointData()->GetScalars();
	vtksys_ios::ostringstream newname;
	newname << (outArray->GetName()?outArray->GetName():"")
    				<< "Slice";
	outArray->SetName(newname.str().c_str());
	// Why not pass the original array?
	if (this->GetInputArrayToProcess(0, inputVector))
	{
		output->GetPointData()->AddArray(
				this->GetInputArrayToProcess(0, inputVector));
	}
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

	// Gradient makes sense only with one input component.  This is not
	// a Jacobian filter.
	if(inputArray->GetNumberOfComponents() != 1)
	{
		vtkErrorMacro(
				"Execute: input has more than one component. "
				"The input to gradient should be a single component image. "
				"Think about it. If you insist on using a color image then "
				"run it though RGBToHSV then ExtractComponents to get the V "
				"components. That's probably what you want anyhow.");
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
	if (Curve != NULL)
		time = std::max(time, Curve->GetMTime());
	if (CutFunction != NULL)
		time = std::max(time, CutFunction->GetMTime());
	return time;
}
