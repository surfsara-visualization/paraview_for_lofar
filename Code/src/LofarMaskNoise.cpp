#include "LofarMaskNoise.h"

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <cmath>
#include <fstream>
#include <vtksys/ios/sstream>

vtkStandardNewMacro(LofarMaskNoise);

/// The function where the actual work is done is LofarIntegrateFrequenciesExecute
/// See Lofar Gaussian blur for more comments on the structure of the code

//----------------------------------------------------------------------------
// Construct an instance of LofarMaskNoise fitler.
LofarMaskNoise::LofarMaskNoise()
{
	// by default process active point scalars
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
			vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
void LofarMaskNoise::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int LofarMaskNoise::RequestInformation(vtkInformation*,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
	// Get input and output pipeline information.
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// Get the input whole extent.
	int extent[6];
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);

	// Store the new whole extent for the output.
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

	// Set the number of point data componets to the number of
	// components in the gradient vector.
	vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_CHAR, 1);

	return 1;
}

//----------------------------------------------------------------------------
// This method computes the input extent necessary to generate the output.
int LofarMaskNoise::RequestUpdateExtent(vtkInformation*,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector)
{
	// Get input and output pipeline information.
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

	// Get the input whole extent.
	int wholeExtent[6];
	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);

	// Get the requested update extent from the output.
	int inUExt[6];
	outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inUExt);
    inUExt[0] = wholeExtent[0];
    inUExt[1] = wholeExtent[1];
    inUExt[2] = wholeExtent[2];
    inUExt[3] = wholeExtent[3];

	// Store the update extent needed from the intput.
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inUExt, 6);

	return 1;
}

//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values
// out of extent.
template <class T>
void LofarMaskNoiseExecute(LofarMaskNoise *self,
		vtkImageData *inData, T *inPtr,
		vtkImageData *outData, char *outPtr,
		int outExt[6], int id)
{
    T *imagePtr = inPtr;
	int idxX, idxY, idxZ;
	int maxX, maxY, maxZ;
	vtkIdType inIncX, inIncY, inIncZ;
	vtkIdType outIncX, outIncY, outIncZ;
	unsigned long count = 0;
	unsigned long target;
	int *inExt = inData->GetExtent();
	int *wholeExtent;
	vtkIdType *inIncs;
	double r[3];

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
	wholeExtent = inData->GetExtent();

    // Move the pointer to the correct starting position.
    imagePtr += (outExt[4]-inExt[4])*inIncs[2];

    // Move the pointer to the correct starting position.
    inPtr = imagePtr + (outExt[0]-inExt[0])*inIncs[0] + (outExt[2]-inExt[2])*inIncs[1];

	// Loop through output pixels
	for (idxZ = wholeExtent[4]; idxZ <= maxZ; idxZ++)
	{
		double sq_average = 0.0;
		double average = 0.0;

		T *inPtr2 = imagePtr;

		// Loop over the entire image to get the distribution

		for (idxY = wholeExtent[2]; !self->AbortExecute && idxY <= wholeExtent[3]; idxY++)
		{
			double sq_average_x = 0.0;
			double average_x = 0.0;
			for (idxX = wholeExtent[0]; idxX <= wholeExtent[1]; idxX++)
			{
				sq_average_x += (*inPtr2)*(*inPtr2);
				average_x += *inPtr2;
				++inPtr2;
			}

			sq_average += sq_average_x;
			average += average_x;
		}

		double n = (wholeExtent[1]-wholeExtent[0])*1.0*(wholeExtent[3]-wholeExtent[2]);
		double mean = average/n;
		double deviation = sqrt(sq_average/n - mean*mean);
		double limit = self->GetStandardDev()*deviation;

		for (idxY = 0; !self->AbortExecute && idxY <= maxY; idxY++)
		{
			if (!id)
			{
				if (!(count%target)) self->UpdateProgress(count/(50.0*target));
				count++;
			}
			for (idxX = 0; idxX <= maxX; idxX++)
			{
				*outPtr = char(*inPtr < limit ? 0 : 1);
				outPtr++;
				inPtr++;
			}
			outPtr += outIncY;
			inPtr += inIncY;
		}
		outPtr += outIncZ;
		inPtr += inIncZ;
        imagePtr += inIncZ;
	}
}

int LofarMaskNoise::RequestData(
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
	newname << (outArray->GetName()?outArray->GetName():"") << "Mask";
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
void LofarMaskNoise::ThreadedRequestData(vtkInformation*,
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
	if(output->GetScalarType() != VTK_CHAR)
	{
		vtkErrorMacro("Execute: output ScalarType is "
				<< output->GetScalarType() << "but must be char.");
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
	char* outPtr = static_cast<char *>(
			output->GetScalarPointerForExtent(outExt));
	switch(inputArray->GetDataType())
	{
	vtkTemplateMacro(
			LofarMaskNoiseExecute(this, input, static_cast<VTK_TT*>(inPtr),
					output, outPtr, outExt, threadId)
	);
	default:
		vtkErrorMacro("Execute: Unknown ScalarType " << input->GetScalarType());
		return;
	}
}
