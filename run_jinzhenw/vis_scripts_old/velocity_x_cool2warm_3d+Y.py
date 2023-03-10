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
    renderView1 = CreateView('RenderView')
    renderView1.ViewSize = [800, 800]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.CenterOfRotation = [256.0, 256.0, 256.0]
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraPosition = [256.0, -1457.1836768696405, 256.0]
    renderView1.CameraFocalPoint = [256.0, 256.0, 256.0]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraViewAngle = 21.849865951742625
    renderView1.CameraFocalDisk = 1.0
    renderView1.CameraParallelScale = 443.4045633326258

    SetActiveView(None)

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.AssignView(0, renderView1)
    layout1.SetSize(800, 800)

    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(renderView1)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'VisItPixieReader'
    nVB_C009_l10n512_S12345T692_z54h5 = VisItPixieReader(registrationName=regName, FileName=inputFilename)
    nVB_C009_l10n512_S12345T692_z54h5.Meshes = ['mesh_512x512x512']
    nVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/velocity_x']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=nVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsvelocity_x'
    native_fieldsvelocity_xLUT = GetColorTransferFunction('native_fieldsvelocity_x')
    native_fieldsvelocity_xLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldsvelocity_xLUT.RGBPoints = [-16080818.828125, 0.0, 0.0, 0.34902, -15161180.631591797, 0.039216, 0.062745, 0.380392, -14241542.435058594, 0.062745, 0.117647, 0.411765, -13321904.23852539, 0.090196, 0.184314, 0.45098, -12402266.041992188, 0.12549, 0.262745, 0.501961, -11482627.845458984, 0.160784, 0.337255, 0.541176, -10562989.648925781, 0.2, 0.396078, 0.568627, -9643351.452392578, 0.239216, 0.454902, 0.6, -8723713.255859375, 0.286275, 0.521569, 0.65098, -7804075.059326172, 0.337255, 0.592157, 0.701961, -6884436.862792969, 0.388235, 0.654902, 0.74902, -5964798.666259766, 0.466667, 0.737255, 0.819608, -5045160.4697265625, 0.572549, 0.819608, 0.878431, -4125522.2731933594, 0.654902, 0.866667, 0.909804, -3205884.0766601562, 0.752941, 0.917647, 0.941176, -2286245.880126953, 0.823529, 0.956863, 0.968627, -1366607.68359375, 0.988235, 0.960784, 0.901961, -1366607.68359375, 0.941176, 0.984314, 0.988235, -778039.2378124986, 0.988235, 0.945098, 0.85098, -189470.79203124903, 0.980392, 0.898039, 0.784314, 472668.70947265625, 0.968627, 0.835294, 0.698039, 1392306.9060058594, 0.94902, 0.733333, 0.588235, 2311945.1025390625, 0.929412, 0.65098, 0.509804, 3231583.2990722656, 0.909804, 0.564706, 0.435294, 4151221.4956054688, 0.878431, 0.458824, 0.352941, 5070859.692138672, 0.839216, 0.388235, 0.286275, 5990497.888671875, 0.760784, 0.294118, 0.211765, 6910136.085205078, 0.701961, 0.211765, 0.168627, 7829774.281738281, 0.65098, 0.156863, 0.129412, 8749412.478271484, 0.6, 0.094118, 0.094118, 9669050.674804688, 0.54902, 0.066667, 0.098039, 10588688.87133789, 0.501961, 0.05098, 0.12549, 11508327.067871094, 0.45098, 0.054902, 0.172549, 12427965.264404297, 0.4, 0.054902, 0.192157, 13347603.4609375, 0.34902, 0.070588, 0.211765]
    native_fieldsvelocity_xLUT.ShowDataHistogram = 1
    native_fieldsvelocity_xLUT.AutomaticDataHistogramComputation = 1
    native_fieldsvelocity_xLUT.ColorSpace = 'Lab'
    native_fieldsvelocity_xLUT.NanColor = [0.25, 0.0, 0.0]
    native_fieldsvelocity_xLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsvelocity_x'
    native_fieldsvelocity_xPWF = GetOpacityTransferFunction('native_fieldsvelocity_x')
    native_fieldsvelocity_xPWF.Points = [-16080818.828125, 0.0, 0.5, 0.0, 13347603.4609375, 1.0, 0.5, 0.0]
    native_fieldsvelocity_xPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Volume'
    resampleToImage1Display.ColorArrayName = ['POINTS', 'native_fields/velocity_x']
    resampleToImage1Display.LookupTable = native_fieldsvelocity_xLUT
    resampleToImage1Display.SelectTCoordArray = 'None'
    resampleToImage1Display.SelectNormalArray = 'None'
    resampleToImage1Display.SelectTangentArray = 'None'
    resampleToImage1Display.OSPRayScaleArray = 'native_fields/velocity_x'
    resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleToImage1Display.SelectOrientationVectors = 'None'
    resampleToImage1Display.ScaleFactor = 51.199948799999994
    resampleToImage1Display.SelectScaleArray = 'None'
    resampleToImage1Display.GlyphType = 'Arrow'
    resampleToImage1Display.GlyphTableIndexArray = 'None'
    resampleToImage1Display.GaussianRadius = 2.5599974399999996
    resampleToImage1Display.SetScaleArray = ['POINTS', 'native_fields/velocity_x']
    resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.OpacityArray = ['POINTS', 'native_fields/velocity_x']
    resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
    resampleToImage1Display.ScalarOpacityUnitDistance = 1.7354386040415883
    resampleToImage1Display.ScalarOpacityFunction = native_fieldsvelocity_xPWF
    resampleToImage1Display.OpacityArrayName = ['POINTS', 'native_fields/velocity_x']
    resampleToImage1Display.SliceFunction = 'Plane'
    resampleToImage1Display.Slice = 255

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage1Display.ScaleTransferFunction.Points = [-27192590.0, 0.0, 0.5, 0.0, 30000606.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage1Display.OpacityTransferFunction.Points = [-27192590.0, 0.0, 0.5, 0.0, 30000606.0, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    resampleToImage1Display.SliceFunction.Origin = [256.0, 256.0, 256.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for native_fieldsvelocity_xLUT in view renderView1
    native_fieldsvelocity_xLUTColorBar = GetScalarBar(native_fieldsvelocity_xLUT, renderView1)
    native_fieldsvelocity_xLUTColorBar.Title = 'native_fields/velocity_x'
    native_fieldsvelocity_xLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldsvelocity_xLUTColorBar.Visibility = 1

    # show color legend
    resampleToImage1Display.SetScalarBarVisibility(renderView1, True)

    SaveScreenshot(outputName, renderView1, ImageResolution=[800, 800], CompressionLevel='0')

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
