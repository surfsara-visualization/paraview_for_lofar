#include "LofarGaussianBlur.h"

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

vtkStandardNewMacro(LofarGaussianBlur);

/// The function where the actual work is done is LofarGaussianBlurExecute
/// See Lofar Gaussian blur for more comments on the structure of the code


//----------------------------------------------------------------------------
// Construct an instance of LofarGaussianBlur filter.
LofarGaussianBlur::LofarGaussianBlur()
{
    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                    vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
void LofarGaussianBlur::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int LofarGaussianBlur::RequestInformation(vtkInformation*,
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

    // Set the number of point data components to one double.
    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);

    return 1;
}

//----------------------------------------------------------------------------
// This method computes the input extent necessary to generate the output.
int LofarGaussianBlur::RequestUpdateExtent(vtkInformation*,
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

    // Store the update extent needed from the input.
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inUExt, 6);

    return 1;
}

//----------------------------------------------------------------------------
// Actual implementation that does the gaussian blur
template <class T>
void LofarGaussianBlurExecute(LofarGaussianBlur *self,
                vtkImageData *inData, T *inPtr,
                vtkImageData *outData, double *outPtr,
                int outExt[6], int id)
{
    int idxZ;
    int maxX, maxY, maxZ;
    vtkIdType inIncX, inIncY, inIncZ;
    vtkIdType outIncX, outIncY, outIncZ;
    unsigned long count = 0;
    unsigned long target;
    int *inExt = inData->GetExtent();
    int *wholeExtent;
    vtkIdType *inIncs;

    // find the region to loop over
    maxX = outExt[1] - outExt[0];
    maxY = outExt[3] - outExt[2];
    maxZ = outExt[5] - outExt[4];
    target = static_cast<unsigned long>((maxZ+1)*(maxY+1)/50.0);
    target++;

    const int N = self->GetKernelSize();
    const int N_over_2 = N/2;
    std::vector<double> kernel(N, 0.0);
    std::vector<double> tmp_img((maxY+1) * (maxX+1), 0.0);
    std::vector<T *>    tmp_rows(N, NULL);

    // Get increments to march through data
    inData->GetContinuousIncrements(outExt, inIncX, inIncY, inIncZ);
    outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

    // get some other info we need
    inIncs = inData->GetIncrements();
    wholeExtent = inData->GetExtent();

    // Move the pointer to the correct starting position.
    inPtr += (outExt[0]-inExt[0])*inIncs[0] +
                    (outExt[2]-inExt[2])*inIncs[1] +
                    (outExt[4]-inExt[4])*inIncs[2];

    {   // The Gaussian kernel
        for (int i=0; i<N; ++i) kernel[i] = pow(M_E, -(i-(N-1)/2)*(i-(N-1)/2)/(2*N*N));
        // Normalize
        double sum = 0.0;
        for (int i=0; i<N; ++i) sum += kernel[i];
        for (int i=0; i<N; ++i) kernel[i] /= sum;
    }

    // Loop through output pixels
    for (idxZ = wholeExtent[4]; idxZ <= maxZ; idxZ++)
    {
        if (!id)
        {
            if (!(count%target)) self->UpdateProgress(count/(50.0*target));
            count++;
        }

        // We split the Gaussian kernel in two parts: in the y-direction and the x-direction
        // This is allowed for a Gaussian kernel and saves calculations

        {   // Copy data to the tmp buffer, y-kernel
            double *tmp_idx = &tmp_img[0];
            for (int y=0; y<=maxY; ++y) {
                // Zero out the values
                for (int x=0; x<=maxX; ++x) tmp_idx[x] = 0.0;

                // Get the rows that are multiplied with the kernel
                for (int yp = std::max(-N_over_2, -y); yp < std::min(maxY-y, N_over_2); ++yp)
                    tmp_rows[N_over_2+yp] = &inPtr[yp * (maxX+1)];

                for (int x=0; x<=maxX; ++x) {
                    // Iterate over the kernel. Special care is needed for the boundary where you cannot use the entire kernel.
                    for (int yp = std::max(y-N_over_2, 0)-y; yp < std::min(maxY, y+N_over_2)-y; ++yp) {
                        // Multiply the kernel with the data
                        *tmp_idx += kernel[N_over_2 + yp] * *tmp_rows[N_over_2+yp];
                        // Go to the next element in the row
                        tmp_rows[N_over_2+yp]++;
                    }

                    // next pixel
                    tmp_idx++;
                    inPtr++;
                }
                inPtr += inIncY;
            }
            inPtr += inIncZ;
        }

        {   // Copy the data to the output buffer, x-kernel
            double *tmp_idx = &tmp_img[0];
            for (int y=0; y<=maxY; ++y) {
                for (int x=0; x<=maxX; ++x) outPtr[x] = 0.0;

                for (int x=0; x<=maxX; ++x) {
                    for (int i=-std::min(x, N_over_2); i<=N_over_2 && x+i <= maxX; ++i)
                        outPtr[i] += kernel[i+N_over_2] * (*tmp_idx);
                    //*outPtr = *tmp_idx;
                    tmp_idx++;
                    outPtr++;
                }
                outPtr += outIncY;
            }
            outPtr += outIncZ;
        }
    }
}

int LofarGaussianBlur::RequestData(
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
            		                << "Blurred";
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
// templated function for the input data type.
void LofarGaussianBlur::ThreadedRequestData(vtkInformation*,
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

    // Gaussian kernel is only implemented for one input component.
    if(inputArray->GetNumberOfComponents() != 1)
    {
        vtkErrorMacro(
                        "Execute: input has more than one component. "
                        "The Gaussian kernel is only implemented for one component.");
        return;
    }

    void* inPtr = inputArray->GetVoidPointer(0);
    double* outPtr = static_cast<double *>(
                    output->GetScalarPointerForExtent(outExt));
    switch(inputArray->GetDataType())
    {
    vtkTemplateMacro(
                    LofarGaussianBlurExecute(this, input, static_cast<VTK_TT*>(inPtr),
                                    output, outPtr, outExt, threadId)
    );
    default:
        vtkErrorMacro("Execute: Unknown ScalarType " << input->GetScalarType());
        return;
    }
}
