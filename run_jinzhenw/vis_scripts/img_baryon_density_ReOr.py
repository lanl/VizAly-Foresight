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
    renderView2.CameraPosition = [171.91689193625075, -1446.9557800825496, 89.05703754169178]
    renderView2.CameraFocalPoint = [256.0, 256.0, 256.0]
    renderView2.CameraViewUp = [0.0, 0.0, 1.0]
    renderView2.CameraViewAngle = 23.637176050044683
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
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/baryon_density']

    # create a new 'Resample To Image'
    resampleToImage2 = ResampleToImage(registrationName='ResampleToImage2', Input=mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage2.SamplingDimensions = [512, 512, 512]
    resampleToImage2.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView2'
    # ----------------------------------------------------------------

    # show data from resampleToImage2
    resampleToImage2Display = Show(resampleToImage2, renderView2, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityLUT = GetColorTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldsbaryon_densityLUT.RGBPoints = [0.0, 0.498039, 0.0, 0.0, 31.3725, 0.6004, 0.0, 0.0, 62.745, 0.702514, 0.000738, 0.000477, 94.11775, 0.773379, 0.095225, 0.061499, 125.49025, 0.843875, 0.189865, 0.12283, 156.86275, 0.891119, 0.294195, 0.203537, 188.23524999999998, 0.937855, 0.397924, 0.283137, 219.60774999999998, 0.963445, 0.476663, 0.316601, 250.9803925, 0.988297, 0.555771, 0.351665, 282.353, 0.990265, 0.646321, 0.436309, 313.7255, 0.992157, 0.735256, 0.519646, 345.098, 0.992157, 0.784468, 0.570827, 376.4705, 0.992249, 0.833218, 0.623483, 407.84325, 0.994218, 0.872587, 0.706159, 439.21575, 0.996186, 0.911419, 0.788189, 470.58825, 0.998155, 0.940946, 0.859054, 500.0, 1.0, 0.968627, 0.92549]
    native_fieldsbaryon_densityLUT.ShowDataHistogram = 1
    native_fieldsbaryon_densityLUT.ColorSpace = 'Lab'
    native_fieldsbaryon_densityLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityPWF = GetOpacityTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityPWF.Points = [0.0, 0.0, 0.5, 0.0, 10.025062561035156, 0.16577540338039398, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]
    native_fieldsbaryon_densityPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage2Display.Representation = 'Volume'
    resampleToImage2Display.ColorArrayName = ['POINTS', 'native_fields/baryon_density']
    resampleToImage2Display.LookupTable = native_fieldsbaryon_densityLUT
    resampleToImage2Display.SelectTCoordArray = 'None'
    resampleToImage2Display.SelectNormalArray = 'None'
    resampleToImage2Display.SelectTangentArray = 'None'
    resampleToImage2Display.OSPRayScaleArray = 'native_fields/baryon_density'
    resampleToImage2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleToImage2Display.SelectOrientationVectors = 'None'
    resampleToImage2Display.ScaleFactor = 51.199948799999994
    resampleToImage2Display.SelectScaleArray = 'None'
    resampleToImage2Display.GlyphType = 'Arrow'
    resampleToImage2Display.GlyphTableIndexArray = 'None'
    resampleToImage2Display.GaussianRadius = 2.5599974399999996
    resampleToImage2Display.SetScaleArray = ['POINTS', 'native_fields/baryon_density']
    resampleToImage2Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleToImage2Display.OpacityArray = ['POINTS', 'native_fields/baryon_density']
    resampleToImage2Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleToImage2Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleToImage2Display.PolarAxes = 'PolarAxesRepresentation'
    resampleToImage2Display.ScalarOpacityUnitDistance = 1.7354386040415883
    resampleToImage2Display.ScalarOpacityFunction = native_fieldsbaryon_densityPWF
    resampleToImage2Display.OpacityArrayName = ['POINTS', 'native_fields/baryon_density']
    resampleToImage2Display.SliceFunction = 'Plane'
    resampleToImage2Display.Slice = 255

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage2Display.ScaleTransferFunction.Points = [-42.67316436767578, 0.0, 0.5, 0.0, 28917.775390625, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage2Display.OpacityTransferFunction.Points = [-42.67316436767578, 0.0, 0.5, 0.0, 28917.775390625, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    resampleToImage2Display.SliceFunction.Origin = [256.0, 256.0, 256.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for native_fieldsbaryon_densityLUT in view renderView2
    native_fieldsbaryon_densityLUTColorBar = GetScalarBar(native_fieldsbaryon_densityLUT, renderView2)
    native_fieldsbaryon_densityLUTColorBar.Title = 'native_fields/baryon_density'
    native_fieldsbaryon_densityLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldsbaryon_densityLUTColorBar.Visibility = 1

    # show color legend
    resampleToImage2Display.SetScalarBarVisibility(renderView2, True)

    SaveScreenshot(outputName, renderView2, ImageResolution=[800, 800], CompressionLevel='0')

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(resampleToImage2)
    # ----------------------------------------------------------------


if __name__ == '__main__':
    print("filename:", sys.argv[1])
    print("registration name:", sys.argv[2])
    print("output filename:", sys.argv[3])

    saveScreenshot(sys.argv[1], sys.argv[2], sys.argv[3])