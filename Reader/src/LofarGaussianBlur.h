#ifndef __LOFARGAUSSIANBLUR_H
#define __LOFARGAUSSIANBLUR_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"

class VTKIMAGINGGENERAL_EXPORT LofarGaussianBlur : public vtkThreadedImageAlgorithm
{
public:
  static LofarGaussianBlur *New();
  vtkTypeMacro(LofarGaussianBlur,vtkThreadedImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetMacro(KernelSize, int);
  vtkGetMacro(KernelSize, int);

protected:
  LofarGaussianBlur();
  ~LofarGaussianBlur() {};

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
  LofarGaussianBlur(const LofarGaussianBlur&);  // Not implemented.
  void operator=(const LofarGaussianBlur&);  // Not implemented.

  int KernelSize;
};

#endif
