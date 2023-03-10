# state file generated using paraview version 5.10.1
import sys

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *

def saveScreenshot(inputFilename, regName, outputName):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # ----------------------------------------------------------------
    # setup views used in the visualization
    # ----------------------------------------------------------------

    # Create a new 'Render View'
    renderView2 = CreateView('RenderView')
    renderView2.ViewSize = [800, 800]
    renderView2.AxesGrid = 'GridAxes3DActor'
    renderView2.CenterOfRotation = [256.0, 256.0, 256.0]
    renderView2.StereoType = 'Crystal Eyes'
    renderView2.CameraPosition = [256.0, -1457.1836768696405, 256.0]
    renderView2.CameraFocalPoint = [256.0, 256.0, 256.0]
    renderView2.CameraViewUp = [0.0, 0.0, 1.0]
    renderView2.CameraViewAngle = 21.849865951742625
    renderView2.CameraFocalDisk = 1.0
    renderView2.CameraParallelScale = 443.4045633326258

    SetActiveView(None)

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.AssignView(0, renderView2)
    layout1.SetSize(800, 800)

    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(renderView2)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'VisItPixieReader'
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5 = VisItPixieReader(registrationName=regName, FileName=inputFilename)
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5.Meshes = ['mesh_512x512x512']
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/temperature']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView2'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView2, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldstemperature'
    native_fieldstemperatureLUT = GetColorTransferFunction('native_fieldstemperature')
    native_fieldstemperatureLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldstemperatureLUT.RGBPoints = [0.0, 0.019608, 0.188235, 0.380392, 6274.499999999999, 0.088504, 0.321107, 0.564937, 12548.999999999998, 0.163399, 0.444983, 0.697501, 18823.55, 0.247059, 0.555709, 0.754095, 25098.05, 0.420684, 0.676432, 0.818685, 31372.55, 0.606459, 0.789773, 0.880277, 37647.049999999996, 0.761476, 0.868512, 0.924567, 43921.549999999996, 0.878047, 0.925721, 0.951942, 50196.078499999996, 0.969089, 0.966474, 0.964937, 56470.600000000006, 0.983852, 0.897578, 0.846828, 62745.1, 0.982468, 0.800692, 0.706113, 69019.6, 0.960323, 0.66782, 0.536332, 75294.1, 0.894579, 0.503806, 0.399769, 81568.65, 0.81707, 0.33218, 0.281046, 87843.15000000001, 0.728489, 0.155017, 0.197386, 94117.65000000001, 0.576932, 0.055363, 0.14925, 100000.0, 0.403922, 0.0, 0.121569]
    native_fieldstemperatureLUT.ShowDataHistogram = 1
    native_fieldstemperatureLUT.AutomaticDataHistogramComputation = 1
    native_fieldstemperatureLUT.ColorSpace = 'Lab'
    native_fieldstemperatureLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldstemperature'
    native_fieldstemperaturePWF = GetOpacityTransferFunction('native_fieldstemperature')
    native_fieldstemperaturePWF.Points = [0.0, 0.0, 0.5, 0.0, 100000.0, 1.0, 0.5, 0.0]
    native_fieldstemperaturePWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Volume'
    resampleToImage1Display.ColorArrayName = ['POINTS', 'native_fields/temperature']
    resampleToImage1Display.LookupTable = native_fieldstemperatureLUT
    resampleToImage1Display.SelectTCoordArray = 'None'
    resampleToImage1Display.SelectNormalArray = 'None'
    resampleToImage1Display.SelectTangentArray = 'None'
    resampleToImage1Display.OSPRayScaleArray = 'native_fields/temperature'
    resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleToImage1Display.SelectOrientationVectors = 'None'
    resampleToImage1Display.ScaleFactor = 51.199948799999994
    resampleToImage1Display.SelectScaleArray = 'None'
    resampleToImage1Display.GlyphType = 'Arrow'
    resampleToImage1Display.GlyphTableIndexArray = 'None'
    resampleToImage1Display.GaussianRadius = 2.5599974399999996
    resampleToImage1Display.SetScaleArray = ['POINTS', 'native_fields/temperature']
    resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.OpacityArray = ['POINTS', 'native_fields/temperature']
    resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
    resampleToImage1Display.ScalarOpacityUnitDistance = 8.957667946113652
    resampleToImage1Display.ScalarOpacityFunction = native_fieldstemperaturePWF
    resampleToImage1Display.OpacityArrayName = ['POINTS', 'native_fields/temperature']
    resampleToImage1Display.SliceFunction = 'Plane'
    resampleToImage1Display.Slice = 49

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage1Display.ScaleTransferFunction.Points = [1194.308349609375, 0.0, 0.5, 0.0, 330162.90625, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage1Display.OpacityTransferFunction.Points = [1194.308349609375, 0.0, 0.5, 0.0, 330162.90625, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    resampleToImage1Display.SliceFunction.Origin = [256.0, 256.0, 256.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for native_fieldstemperatureLUT in view renderView2
    native_fieldstemperatureLUTColorBar = GetScalarBar(native_fieldstemperatureLUT, renderView2)
    native_fieldstemperatureLUTColorBar.Title = 'native_fields/temperature'
    native_fieldstemperatureLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldstemperatureLUTColorBar.Visibility = 1

    # show color legend
    resampleToImage1Display.SetScalarBarVisibility(renderView2, True)

    SaveScreenshot(outputName, renderView2, ImageResolution=[800, 800], CompressionLevel='0')

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(resampleToImage1)
    # ----------------------------------------------------------------


if __name__ == '__main__':
    print("filename:", sys.argv[1])
    print("registration name:", sys.argv[2])
    print("output filename:", sys.argv[3])

    saveScreenshot(sys.argv[1], sys.argv[2], sys.argv[3])
