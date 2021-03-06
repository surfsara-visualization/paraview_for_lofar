#include "LofarIntegrateFrequencies.h"

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>
#include <fstream>
#include <vtksys/ios/sstream>

vtkStandardNewMacro(LofarIntegrateFrequencies);

/// The function where the actual work is done is LofarIntegrateFrequenciesExecute
/// See Lofar Gaussian blur for more comments on the structure of the code

//----------------------------------------------------------------------------
// Construct an instance of LofarIntegrateFrequencies fitler.
LofarIntegrateFrequencies::LofarIntegrateFrequencies()
{
    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
void LofarIntegrateFrequencies::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int LofarIntegrateFrequencies::RequestInformation(vtkInformation*,
                vtkInformationVector** inputVector,
                vtkInformationVector* outputVector)
{
    // Get input and output pipeline information.
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

    // Get the input whole extent.
    int extent[6];
    inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
    extent[4] = 0;
    extent[5] = 0;

    // Store the new whole extent for the output.
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

    // Set the number of point data components to one double.
    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);

    return 1;
}

//----------------------------------------------------------------------------
// This method computes the input extent necessary to generate the output.
int LofarIntegrateFrequencies::RequestUpdateExtent(vtkInformation*,
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
    inUExt[4] = wholeExtent[4];
    inUExt[5] = wholeExtent[5];

    // Store the update extent needed from the intput.
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inUExt, 6);

    return 1;
}

template <class T>
void LofarIntegrateFrequenciesExecute(LofarIntegrateFrequencies *self,
                vtkImageData *inData, T *inPtr,
                vtkImageData *outData, double *outPtr,
                int outExt[6], int id)
{
    int idxX, idxY, idxZ;
    int maxX, maxY;
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

    target = static_cast<unsigned long>((wholeExtent[5]-wholeExtent[4]+1)*(maxY+1)/50.0);
    target++;

    // Move the pointer to the correct starting position.
    inPtr += (outExt[0]-inExt[0])*inIncs[0] +
                    (outExt[2]-inExt[2])*inIncs[1] +
                    (outExt[4]-inExt[4])*inIncs[2];

    // Zero the output image
    double *outPtr2 = outPtr;
    for (idxY = 0; idxY <= maxY; idxY++) {
        for (idxX = 0; idxX <= maxX; idxX++) {
            *outPtr2 = 0.0;
            outPtr2++;
        }
    }

    // Loop through output pixels
    for (idxZ = wholeExtent[4]; idxZ <= wholeExtent[5]; idxZ++)
    {
        outPtr2 = outPtr;
        for (idxY = 0; !self->AbortExecute && idxY <= maxY; idxY++)
        {
            if (!id)
            {
                if (!(count%target)) {
                    self->UpdateProgress(count/(50.0*target));
                }
                count++;
            }
            for (idxX = 0; idxX <= maxX; idxX++)
            {
                // Just sum up every frequency ...
                *outPtr2 += *inPtr;
                outPtr2++;
                inPtr++;
            }
            outPtr2 += outIncY;
            inPtr += inIncY;
        }
        inPtr += inIncZ;
    }
}

int LofarIntegrateFrequencies::RequestData(
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
    newname << (outArray->GetName()?outArray->GetName():"") << "Accumulated";
    outArray->SetName(newname.str().c_str());

    return 1;
}

//----------------------------------------------------------------------------
// This method contains a switch statement that calls the correct
// templated function for the input data type.  This method does handle
// boundary conditions.
void LofarIntegrateFrequencies::ThreadedRequestData(vtkInformation*,
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

    // Integration only works on scalars.
    if(inputArray->GetNumberOfComponents() != 1)
    {
        vtkErrorMacro("Execute: input has more than one component. ");
        return;
    }

    void* inPtr = inputArray->GetVoidPointer(0);
    double* outPtr = static_cast<double *>(
                    output->GetScalarPointerForExtent(outExt));
    switch(inputArray->GetDataType())
    {
    vtkTemplateMacro(
                    LofarIntegrateFrequenciesExecute(this, input, static_cast<VTK_TT*>(inPtr),
                                    output, outPtr, outExt, threadId)
    );
    default:
        vtkErrorMacro("Execute: Unknown ScalarType " << input->GetScalarType());
        return;
    }
}
