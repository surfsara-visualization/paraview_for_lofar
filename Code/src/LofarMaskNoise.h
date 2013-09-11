#ifndef __LofarMaskNoise_H
#define __LofarMaskNoise_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"

class VTKIMAGINGGENERAL_EXPORT LofarMaskNoise : public vtkThreadedImageAlgorithm
{
public:
	static LofarMaskNoise *New();
	vtkTypeMacro(LofarMaskNoise,vtkThreadedImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	vtkSetMacro(StandardDev, double);
	vtkGetMacro(StandardDev, double);

protected:
	LofarMaskNoise();
	~LofarMaskNoise() {};

	virtual int RequestInformation (vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);
	virtual int RequestUpdateExtent(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);
	virtual int RequestData(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*);

	void ThreadedRequestData(vtkInformation*,
			vtkInformationVector**,
			vtkInformationVector*,
			vtkImageData*** inData,
			vtkImageData** outData,
			int outExt[6],
			int threadId);
private:
	LofarMaskNoise(const LofarMaskNoise&);  // Not implemented.
	void operator=(const LofarMaskNoise&);  // Not implemented.

	double StandardDev;
};

#endif
