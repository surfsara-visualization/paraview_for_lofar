#ifndef __LofarApplyMask_H
#define __LofarApplyMask_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"

#include <iostream>

class VTKIMAGINGGENERAL_EXPORT LofarApplyMask : public vtkThreadedImageAlgorithm
{
public:
	static LofarApplyMask *New();
	vtkTypeMacro(LofarApplyMask,vtkThreadedImageAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

protected:
	LofarApplyMask();
	~LofarApplyMask() {};

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
	LofarApplyMask(const LofarApplyMask&);  // Not implemented.
	void operator=(const LofarApplyMask&);  // Not implemented.
};

#endif
