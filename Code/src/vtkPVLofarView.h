/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile$

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPVLofarView
// .SECTION Description
// vtkPVLofarView demonstrates how to create custom render-view
// subclasses that use a image-processing render pass for processing the image
// before rendering it on the screen.

#ifndef __vtkPVLofarView_h
#define __vtkPVLofarView_h

#include "vtkPVRenderView.h"
#include "vtkTextActor.h"
#include "vtkObjectFactory.h"
#include "vtkPVSynchronizedRenderer.h"


#include "vtk3DWidgetRepresentation.h"
#include "vtkAbstractWidget.h"
#include "vtkCamera.h"
#include "vtkCameraManipulator.h"
#include "vtkCommand.h"
#include "vtkHandleRepresentation.h"
#include "vtkHandleWidget.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkPVGenericRenderWindowInteractor.h"
#include "vtkPVInteractorStyle.h"
#include "vtkPVSynchronizedRenderWindows.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkPointSource.h"
#include "vtkRenderViewBase.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkWeakPointer.h"

class VTK_EXPORT vtkPVLofarView : public vtkPVRenderView
{
public:
  static vtkPVLofarView* New();
  vtkTypeMacro(vtkPVLofarView, vtkPVRenderView);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Initialize the view with an identifier. Unless noted otherwise, tshis method
  // must be called before calling any other methods on this class.
  // @CallOnAllProcessess
    virtual void Initialize(unsigned int id);
    void SetScalarValue(double value);
    const char * GetXlabel(void);
    void SetXLabel(const char* text);


//BTX
protected:
    vtkPVLofarView();
    ~vtkPVLofarView();

    const char* XAxisLabel;
    char* YAxisLabel;
    char* ZAxisLabel;
    char* ScalarLabel;
    int LabelFontSize;
    double ScalarValue;

    double TextValues[4];
    double NaturalCoordinates[3];
    double camPos[3];
    double camCenter[2];

    vtkNew<vtkTextActor> infoLabel;

    void UpdateDataEventCallBack(vtkObject* src, unsigned long event, void* data);

private:
  vtkPVLofarView(const vtkPVLofarView&); // Not implemented
  void operator=(const vtkPVLofarView&); // Not implemented
//ETX
};

#endif
