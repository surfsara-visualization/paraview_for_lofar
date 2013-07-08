import os.path

try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

#Connect('localhost')

color_lookup_table = \
  GetLookupTableForArray( "intensity", 1, 
    RGBPoints=[0.0,   0.23, 0.299, 0.754,
               0.003, 0.706, 0.016, 0.15],
    #LockScalarRange=1,
    ColorSpace='Diverging'
  )

sources = []
def load_fits(filename):
  obj_name = os.path.splitext(os.path.basename(filename))[0]
  result = FITSreader( FileName=filename, registrationName=obj_name )
  sources.append(result)
  Show()
  return result

data_directory = '/home/nicok/Projects/lofar/Data'
#for x in range(1):
#  for y in range(1):
#    filename = data_directory+'/offset_'+str(x)+'_'+str(y)+'.fits'
#    load_fits(filename)

#load_fits(data_directory+"/ngc4565_30arcsec.fits")
load_fits(data_directory+"/p-rmcube-clean.fits")


def change_representation_to_slice(representation, slice_nr = 0):
    representation.Representation = 'Slice'
    representation.ScalarOpacityFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    representation.ColorArrayName = 'intensity'
    representation.Slice = slice_nr
    representation.LookupTable = color_lookup_table
RenderView1 = GetRenderView()
RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
RenderView1.CameraParallelProjection = 1
RenderView1.ResetCamera()
Render()


for src in sources:
    src_representation = GetDisplayProperties(src)
    change_representation_to_slice(src_representation, 130)
#    Hide()

    slice = Sliceimagedata(Input=src)
    slice_representation = GetDisplayProperties(slice)
    #change_representation_to_slice(slice_representation)
    slice.Slicecurve.Point1 = [3.313, 0.450, 0.003]
    slice.Slicecurve.Point2 = [3.307, 0.457, 0.003]
    Show()

RenderView2 = CreateRenderView()
SetActiveView(RenderView2)
slice_representation = GetDisplayProperties(slice)
#slice_representation.ColorArrayName = 'intensitySlice'
Show()
RenderView2.CameraParallelProjection = 1
RenderView2.ResetCamera()

AnimationScene1 = GetAnimationScene()
AnimationScene1.ViewModules = [ RenderView1, RenderView2 ]
Render()

a1_intensitySlice_PVLookupTable = GetLookupTableForArray( "intensitySlice", 1, RGBPoints=[-0.0009598891483619809, 0.23, 0.299, 0.754, 0.1074233204126358, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_intensitySlice_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

slice_representation.ScalarOpacityFunction = a1_intensitySlice_PiecewiseFunction
slice_representation.ColorArrayName = 'intensitySlice'
slice_representation.LookupTable = a1_intensitySlice_PVLookupTable

Render()
