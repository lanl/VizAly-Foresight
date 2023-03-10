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
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/velocity_x']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView2'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView2, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsvelocity_x'
    native_fieldsvelocity_xLUT = GetColorTransferFunction('native_fieldsvelocity_x')
    native_fieldsvelocity_xLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldsvelocity_xLUT.RGBPoints = [-17007040.7578125, 0.0, 0.0, 0.34902, -15970781.320800781, 0.039216, 0.062745, 0.380392, -14934521.883789062, 0.062745, 0.117647, 0.411765, -13898262.446777344, 0.090196, 0.184314, 0.45098, -12862003.009765625, 0.12549, 0.262745, 0.501961, -11825743.572753906, 0.160784, 0.337255, 0.541176, -10789484.135742188, 0.2, 0.396078, 0.568627, -9753224.698730469, 0.239216, 0.454902, 0.6, -8716965.26171875, 0.286275, 0.521569, 0.65098, -7680705.824707031, 0.337255, 0.592157, 0.701961, -6644446.3876953125, 0.388235, 0.654902, 0.74902, -5608186.950683594, 0.466667, 0.737255, 0.819608, -4571927.513671875, 0.572549, 0.819608, 0.878431, -3535668.0766601562, 0.654902, 0.866667, 0.909804, -2499408.6396484375, 0.752941, 0.917647, 0.941176, -1463149.2026367188, 0.823529, 0.956863, 0.968627, -426889.765625, 0.988235, 0.960784, 0.901961, -426889.765625, 0.941176, 0.984314, 0.988235, 236316.2740624994, 0.988235, 0.945098, 0.85098, 899522.3137500025, 0.980392, 0.898039, 0.784314, 1645629.1083984375, 0.968627, 0.835294, 0.698039, 2681888.5454101562, 0.94902, 0.733333, 0.588235, 3718147.982421875, 0.929412, 0.65098, 0.509804, 4754407.419433594, 0.909804, 0.564706, 0.435294, 5790666.8564453125, 0.878431, 0.458824, 0.352941, 6826926.293457031, 0.839216, 0.388235, 0.286275, 7863185.73046875, 0.760784, 0.294118, 0.211765, 8899445.167480469, 0.701961, 0.211765, 0.168627, 9935704.604492188, 0.65098, 0.156863, 0.129412, 10971964.041503906, 0.6, 0.094118, 0.094118, 12008223.478515625, 0.54902, 0.066667, 0.098039, 13044482.915527344, 0.501961, 0.05098, 0.12549, 14080742.352539062, 0.45098, 0.054902, 0.172549, 15117001.789550781, 0.4, 0.054902, 0.192157, 16153261.2265625, 0.34902, 0.070588, 0.211765]
    native_fieldsvelocity_xLUT.ShowDataHistogram = 1
    native_fieldsvelocity_xLUT.AutomaticDataHistogramComputation = 1
    native_fieldsvelocity_xLUT.ColorSpace = 'Lab'
    native_fieldsvelocity_xLUT.NanColor = [0.25, 0.0, 0.0]
    native_fieldsvelocity_xLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsvelocity_x'
    native_fieldsvelocity_xPWF = GetOpacityTransferFunction('native_fieldsvelocity_x')
    native_fieldsvelocity_xPWF.Points = [-17007040.7578125, 0.0, 0.5, 0.0, -16882378.0, 0.0, 0.5, 0.0, 16153261.2265625, 1.0, 0.5, 0.0]
    native_fieldsvelocity_xPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Volume'
    resampleToImage1Display.ColorArrayName = ['POINTS', 'native_fields/velocity_x']
    resampleToImage1Display.LookupTable = native_fieldsvelocity_xLUT
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
    resampleToImage1Display.ScalarOpacityUnitDistance = 1.7354386040415883
    resampleToImage1Display.ScalarOpacityFunction = native_fieldsvelocity_xPWF
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

    # get color legend/bar for native_fieldsvelocity_xLUT in view renderView2
    native_fieldsvelocity_xLUTColorBar = GetScalarBar(native_fieldsvelocity_xLUT, renderView2)
    native_fieldsvelocity_xLUTColorBar.Title = 'native_fields/velocity_x'
    native_fieldsvelocity_xLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldsvelocity_xLUTColorBar.Visibility = 1

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