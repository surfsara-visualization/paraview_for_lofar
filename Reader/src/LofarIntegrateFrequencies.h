#ifndef __LofarIntegrateFrequencies_H
#define __LofarIntegrateFrequencies_H

#include "vtkImagingGeneralModule.h" // For export macro
#include "vtkThreadedImageAlgorithm.h"

class VTKIMAGINGGENERAL_EXPORT LofarIntegrateFrequencies : public vtkThreadedImageAlgorithm
{
public:
  static LofarIntegrateFrequencies *New();
  vtkTypeMacro(LofarIntegrateFrequencies,vtkThreadedImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  LofarIntegrateFrequencies();
  ~LofarIntegrateFrequencies() {};

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
  LofarIntegrateFrequencies(const LofarIntegrateFrequencies&);  // Not implemented.
  void operator=(const LofarIntegrateFrequencies&);  // Not implemented.
};

#endif
