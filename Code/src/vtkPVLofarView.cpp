/*=========================================================================

LofarView Extension

=========================================================================*/
#include "vtkPVLofarView.h"
#include "vtkObjectFactory.h"
#include "vtk3DWidgetRepresentation.h"
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
#include "vtkRenderViewBase.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkWeakPointer.h"
#include "vtkStringArray.h"
#include <QMouseEvent>
#include <typeinfo>
#include <vtkRenderWindowInteractor.h>
#include <sstream>
#include <vtkImageData.h>
#include <vtkImageSliceRepresentation.h>
#include "vtkFieldData.h"
#include "vtkPointData.h"
 
//#define LOG(msg)
//template <class T>
//void suppress_unused_warning(const T &) {}
#define LOG(msg) vtkWarningMacro(<< __PRETTY_FUNCTION__ << ", l" << __LINE__ << ": " << msg);
//#define LOG(msg) {output_file << __PRETTY_FUNCTION__ << ", l" << __LINE__ << ": " << msg << std::endl;}

vtkStandardNewMacro(vtkPVLofarView);
//----------------------------------------------------------------------------
vtkPVLofarView::vtkPVLofarView()
{
    for(int i=0; i < 3; ++i)
    {
        this->NaturalCoordinates[i] = 0.0;
        this->TextValues[i] = 0.0;
    }
    // Scalar init
    this->TextValues[3] = 0.0;
    
    //this->AddObserver(vtkCommand::UpdateDataEvent, this, &vtkPVLofarView::UpdateDataEventCallBack);
    // Add mouse listen event to the render window interactor
    this->GetInteractor()->AddObserver(vtkCommand::MouseMoveEvent, this, &vtkPVLofarView::UpdateDataEventCallBack);
    
    //add the Coordinate X label actor to the renderer
    this->GetNonCompositedRenderer()->AddActor(this->infoLabel.GetPointer());
    
    //Set properties of  infoLabel
    this->infoLabel->SetInput("No Data");
    this->infoLabel->SetPosition(20,20);
    this->infoLabel->GetTextProperty()->SetFontSize(20);
    //this->SetScalarValue(20.0);
        
    //set view modes
    this->GetActiveCamera()->ParallelProjectionOn();
    this->GetActiveCamera()->GetPosition(camPos);
    this->SetInteractionMode(INTERACTION_MODE_2D);
    this->SetCenterAxesVisibility(false);
    this->SetOrientationAxesVisibility(false);


}

//----------------------------------------------------------------------------
vtkPVLofarView::~vtkPVLofarView()
{
    
}

//----------------------------------------------------------------------------
void vtkPVLofarView::Initialize(unsigned int id)
{
    this->Superclass::Initialize(id);
      
    this->GetActiveCamera()->ParallelProjectionOn();
    this->SetInteractionMode(INTERACTION_MODE_2D);
    this->SetCenterAxesVisibility(false);
    this->SetOrientationAxesVisibility(false);
 
}

//----------------------------------------------------------------------------
void vtkPVLofarView::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}


void vtkPVLofarView::SetScalarValue(double value)
{
  this->ScalarValue = value;
}

void vtkPVLofarView::SetXLabel(const char* text)
{
  this->XAxisLabel = text;
}

const char * vtkPVLofarView::GetXlabel(void)
{
  return(this->XAxisLabel);
}

void vtkPVLofarView::UpdateDataEventCallBack(vtkObject* src, unsigned long event, void* data)
{
    vtkRenderWindowInteractor * iren = (vtkRenderWindowInteractor * ) src;
    int pos[2];
    std::stringstream text; //Text is for our label. During developing doubles as debug print statement
    vtkImageSliceRepresentation* target;
    
    this->SetInteractionMode(INTERACTION_MODE_2D);        //If i dont set it here again to 2D, it seems to reset
    
    //Get mouse even pos
    iren->GetEventPosition(pos);
    text << pos[0] << "," << pos[1] << " ";
    
    this->GetActiveCamera()->GetPosition(camPos);
    this->GetActiveCamera()->GetWindowCenter(camCenter );
    //double x, y, z;
    //x = pos[0];
    //y=pos[1];
    //z=0.5;
    //this->GetRenderer ()-> ViewToWorld(x, y, z );	
    //camPos[0]=x;
    //camPos[1]=y;
    //camPos[2]=z;
    //text << " camP: " << camPos[0] << "," << camPos[1] << "," << camPos[2] <<" ";
    //text << " camC: " << camCenter[0] << "," << camCenter[1] << " ";

    for (int i = 0; i <this->GetNumberOfRepresentations();  i++)
    {
        if(target = vtkImageSliceRepresentation::SafeDownCast(this->GetRepresentation(i)))
        {
            
        
            vtkFieldData* fieldData = target->GetRenderedDataObject(0)->GetFieldData();
            vtkStringArray* titleX =  vtkStringArray::SafeDownCast(fieldData->GetAbstractArray("AxisTitleForX"));
            //text << titleX->GetValue(0).c_str() << " ";
            if (vtkImageData*  image = vtkImageData::SafeDownCast(target->GetRenderedDataObject(0))) 
            { 

                int extent[6];
                image->GetExtent (extent);
                //text << extent[0] <<"," <<extent[1] << ",";
                if (vtkPointData *pd = image->GetPointData())
                {
                    if(vtkDataArray * data =  pd->GetScalars("intensity"))
                    {
                       
                        //text << data->GetElementComponentSize() << ":";
                        double result[4];
                        int number = data->GetNumberOfTuples();
                        int index = pos[0]*512*149+pos[1]*149+100;
                        //text << " " << index <<  " ";
                        if (index < number && pos[0] >= 0 && pos[0] <=511 && pos[1] >=0 && pos[1] <=511 )
                        {
                            data->GetTuple(index, result);
                            //text << "tuples:" << number << ":" << result[0] <<","  << result[1] <<","  << result[2] <<","  << result[3];
                        }
                    }
                }
                
            }

        }
    }
    
    //write contents of text to infoLabel
    this->infoLabel->SetInput(text.str().c_str());
   
}



