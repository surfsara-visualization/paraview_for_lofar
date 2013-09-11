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

#include <ctime> 
#include <fftw3.h>
#include <boost/scoped_array.hpp>

vtkStandardNewMacro(LofarGaussianBlur);

//----------------------------------------------------------------------------
// Construct an instance of LofarGaussianBlur fitler.
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

	// Set the number of point data componets to the number of
	// components in the gradient vector.
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

	// Store the update extent needed from the intput.
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inUExt, 6);

	return 1;
}

//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values
// out of extent.
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
	inPtr += (outExt[0]-inExt[0])*inIncs[0] +
			(outExt[2]-inExt[2])*inIncs[1] +
			(outExt[4]-inExt[4])*inIncs[2];

	const int N = self->GetKernelSize();
	const int N_over_2 = N/2;
	boost::scoped_array<double> kernel(new double[N]);
	{
		double sum = 0.0;
		for (int i=0; i<N; ++i) kernel[i] = pow(M_E, -(i-(N-1)/2)*(i-(N-1)/2)/(2*N*N));

		for (int i=0; i<N; ++i) sum += kernel[i];
		for (int i=0; i<N; ++i) kernel[i] /= sum;
	}

	boost::scoped_array<double> tmp_img(new double[(maxY+1) * (maxX+1)]);

	//    std::clock_t duration = 0;
	//    std::clock_t total_start = clock();


	// Loop through output pixels
	for (idxZ = wholeExtent[4]; idxZ <= maxZ; idxZ++)
	{
		if (!id)
		{
			if (!(count%target)) self->UpdateProgress(count/(50.0*target));
			count++;
		}

		{   // Copy data to the tmp buffer, y-kernel
			boost::scoped_array<T *> tmp_rows(new T *[N]);
			double *tmp_idx = tmp_img.get();
			for (int y=0; y<=maxY; ++y) {
				// Zero out the values
				for (int x=0; x<=maxX; ++x) tmp_idx[x] = 0.0;

				// Get the rows that are multiplied with the kernel
				for (int yp = std::max(-N_over_2, -y); yp < std::min(maxY-y, N_over_2); ++yp)
					tmp_rows[N_over_2+yp] = &inPtr[yp * (maxX+1)];

				for (int x=0; x<=maxX; ++x) {
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

		//        std::clock_t start = clock();
		{   // Copy the data to the output buffer, x-kernel
			double *tmp_idx = tmp_img.get();
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
		//        std::clock_t stop = clock();
		//        duration += stop - start;
	}
	//    std::clock_t total_stop = clock();
	//    std::cout << "X     duration : " << duration/(double)CLOCKS_PER_SEC << "s (was 1.8)"  << std::endl;
	//    std::cout << "Total duration : " << (total_stop-total_start)/(double)CLOCKS_PER_SEC << "s (was 3.88)"  << std::endl;
}
//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values
// out of extent.
template <class T>
void LofarGaussianBlurExecuteFFTW(LofarGaussianBlur *self,
		vtkImageData *inData, T *inPtr,
		vtkImageData *outData, double *outPtr,
		int outExt[6], int /*id*/)
{
	int maxX, maxY, maxZ;
	vtkIdType inIncX, inIncY, inIncZ;
	vtkIdType outIncX, outIncY, outIncZ;
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
	inPtr += (outExt[0]-inExt[0])*inIncs[0] +
			(outExt[2]-inExt[2])*inIncs[1] +
			(outExt[4]-inExt[4])*inIncs[2];

	double *spatial_slice = (double *)fftw_malloc((2*maxY * 2*maxX) * sizeof(double));

	fftw_complex *freq_slice = (fftw_complex *)fftw_malloc((2*maxY * 2*maxX) * sizeof(fftw_complex));
	fftw_complex *freq_gauss_kernel = (fftw_complex *)fftw_malloc((2*maxY * 2*maxX) * sizeof(fftw_complex));

	{
		double *idx = spatial_slice;
		for (int y=0; y<2*maxY; ++y) {
			for (int x=0; x<2*maxX; ++x) {
				*idx = 0.0;
				++idx;
			}
		}
	}

	{
		fftw_complex *idx = freq_gauss_kernel;
		for (int y=0; y<2*maxY; ++y) {
			for (int x=0; x<2*maxX; ++x) {
				idx[0][0] = 0;
				idx[0][1] = 0;
				++idx;
			}
		}
		freq_gauss_kernel[0][0] = 1.0;
		freq_gauss_kernel[1][0] = 1.0;
		freq_gauss_kernel[2][0] = 1.0;
		freq_gauss_kernel[3][0] = 1.0;
	}

	// Plan forward image transform
	fftw_plan forward_plan = fftw_plan_dft_r2c_2d(2*maxX, 2*maxY, spatial_slice, freq_slice, FFTW_MEASURE);
	fftw_plan inverse_plan = fftw_plan_dft_c2r_2d(2*maxX, 2*maxY, freq_slice, spatial_slice, FFTW_MEASURE);

	for (int idxZ = wholeExtent[4]; !self->AbortExecute && idxZ <= maxZ; idxZ++)
	{
		{   // Copy data to the fftw buffer
			double *input = spatial_slice;
			for (int y=0; y<maxY; ++y) {
				for (int x=0; x<maxX; ++x) {
					*input = *inPtr;
					input++;
					inPtr++;
				}
				inPtr += inIncY;
			}
			inPtr += inIncZ;
		}

		fftw_execute(forward_plan);

		if (false)
		{   // Multiply (convolution)
			fftw_complex *freq_img_ptr = freq_slice;
			fftw_complex *freq_gauss_ptr = freq_gauss_kernel;
			for (int y=0; y<2*maxY; ++y) {
				for (int x=0; x<2*maxX; ++x) {
					double r = freq_img_ptr[0][0];
					double i = freq_gauss_ptr[0][1];
					freq_slice[0][0] = (r * freq_gauss_ptr[0][0] - i * freq_gauss_ptr[0][1]);
					freq_slice[0][1] = (r * freq_gauss_ptr[0][1] + i * freq_gauss_ptr[0][0]);
					freq_img_ptr ++;
					freq_gauss_ptr ++;
				}
			}
		}

		fftw_execute(inverse_plan);

		{   // Copy the data to the output buffer
			double *input = spatial_slice;
			for (int y=0; y<maxY; ++y) {
				for (int x=0; x<maxX; ++x) {
					*outPtr = 42;//*input;
					input++;
					outPtr++;
				}
				input += maxX;
				outPtr += outIncY;
			}
			outPtr += outIncZ;
		}
	}

	fftw_free(spatial_slice);
	fftw_free(freq_slice);
	fftw_free(freq_gauss_kernel);

	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(inverse_plan);
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
// templated function for the input data type.  This method does handle
// boundary conditions.
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
			LofarGaussianBlurExecute(this, input, static_cast<VTK_TT*>(inPtr),
					output, outPtr, outExt, threadId)
	);
	default:
		vtkErrorMacro("Execute: Unknown ScalarType " << input->GetScalarType());
		return;
	}
}
