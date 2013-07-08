try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

def change_to_slice(rep, range_name):
    color_lookup_table = \
        GetLookupTableForArray( range_name, 1, 
                                RGBPoints=[0.0,   0.23, 0.299, 0.754,
                                           0.001, 0.706, 0.016, 0.15],
                                LockScalarRange=1, ColorSpace='Diverging')
    rep = GetDisplayProperties(src)
    rep.Representation = 'Slice'
    rep.ScalarOpacityFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    rep.ColorArrayName = range_name
    rep.Slice = 20
    rep.LookupTable = color_lookup_table


fits_data = FITSreader( FileName='/home/nicok/Projects/lofar/Data/p-rmcube-clean.fits' )
#fits_data = FITSreader( FileName='/home/nicok/Projects/lofar/Data/ngc4565_30arcsec.fits' )
DataRepresentation1 = Show()

Gaussianblur = Gaussianblur()
DataRepresentation2 = Show()

Masknoise = Masknoise()
DataRepresentation3 = Show()

Applymask = Applymask()
Applymask.SelectallTrackswhere = ['POINTS', 'intensity']
Applymask.withrestriction = ['POINTS', 'intensityBlurredMask']
DataRepresentation4 = Show()

IntegrateZ2 = IntegrateZ()
DataRepresentation5 = Show()


a1_intensity_PVLookupTable = GetLookupTableForArray( "intensity", 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.1362226903438568, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )
a1_intensity_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
a1_intensity_PVLookupTable.ScalarOpacityFunction = a1_intensity_PiecewiseFunction

DataRepresentation1.Representation = 'Slice'
DataRepresentation1.ScalarOpacityFunction = a1_intensity_PiecewiseFunction
DataRepresentation1.ColorArrayName = 'intensity'
DataRepresentation1.LookupTable = a1_intensity_PVLookupTable



a1_intensityAccumulated_PVLookupTable = GetLookupTableForArray( "intensityBlurred", 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.1362226903438568, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_intensityAccumulated_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = ''
DataRepresentation2.ScalarOpacityFunction = a1_intensityAccumulated_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'intensityBlurred'
DataRepresentation2.ScalarOpacityUnitDistance = 0.0005393250329026173
DataRepresentation2.LookupTable = a1_intensityAccumulated_PVLookupTable
DataRepresentation2.Representation = 'Slice'
DataRepresentation2.ScaleFactor = 0.027663468644099055

a1_intensityAccumulated_PVLookupTable.ScalarOpacityFunction = a1_intensityAccumulated_PiecewiseFunction

my_representation1 = GetDisplayProperties(Masknoise)
a1_intensityAccumulatedAccumulated_PVLookupTable = GetLookupTableForArray( "intensityBlurred", 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 1.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )
a1_intensityAccumulatedAccumulated_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

DataRepresentation3.Slice = 50
DataRepresentation3.ScalarOpacityFunction = a1_intensityAccumulatedAccumulated_PiecewiseFunction
DataRepresentation3.ColorArrayName = 'intensityBlurredMask'
DataRepresentation3.Visibility = 1
DataRepresentation3.LookupTable = a1_intensityAccumulatedAccumulated_PVLookupTable
DataRepresentation3.Representation = 'Slice'




DataRepresentation1.Slice = 50
DataRepresentation2.Slice = 50
DataRepresentation3.Slice = 50


DataRepresentation1.Visibility = 0
DataRepresentation2.Visibility = 0
DataRepresentation3.Visibility = 1

RenderView = GetRenderView()
RenderView.ResetCamera()
Render()
