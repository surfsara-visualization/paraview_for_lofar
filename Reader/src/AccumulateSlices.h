#ifndef __ACCUMULATESLICES_H
#define __ACCUMULATESLICES_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"

class VTKIMAGINGGENERAL_EXPORT AccumulateSlices : public vtkThreadedImageAlgorithm
{
public:
  static AccumulateSlices *New();
  vtkTypeMacro(AccumulateSlices,vtkThreadedImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  AccumulateSlices();
  ~AccumulateSlices() {};

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
  AccumulateSlices(const AccumulateSlices&);  // Not implemented.
  void operator=(const AccumulateSlices&);  // Not implemented.
};

#endif
