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
    mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/dark_matter_density']

    # create a new 'Resample To Image'
    resampleToImage2 = ResampleToImage(registrationName='ResampleToImage2', Input=mgard_rel__1e1__NVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage2.SamplingDimensions = [512, 512, 512]
    resampleToImage2.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView2'
    # ----------------------------------------------------------------

    # show data from resampleToImage2
    resampleToImage2Display = Show(resampleToImage2, renderView2, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsdark_matter_density'
    native_fieldsdark_matter_densityLUT = GetColorTransferFunction('native_fieldsdark_matter_density')
    native_fieldsdark_matter_densityLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldsdark_matter_densityLUT.RGBPoints = [3.0, 0.175119, 0.0438468, 1.0, 15.409423999999998, 0.22383, 0.159771, 0.94557, 27.8189465, 0.27254, 0.233611, 0.891216, 40.2283705, 0.321251, 0.296526, 0.836857, 52.637794500000005, 0.369962, 0.354296, 0.782359, 65.04721849999999, 0.418672, 0.409139, 0.72754, 77.456741, 0.467383, 0.462152, 0.672148, 89.866165, 0.51609, 0.51396, 0.615825, 102.27559097000001, 0.572863, 0.55452, 0.559172, 114.685013, 0.630269, 0.593822, 0.517729, 127.0945355, 0.689588, 0.624668, 0.47446, 139.50395949999998, 0.745394, 0.656113, 0.428638, 151.9133835, 0.798624, 0.688104, 0.379105, 164.3228075, 0.849926, 0.720593, 0.323834, 176.73233000000002, 0.899765, 0.753543, 0.258657, 189.141754, 0.948487, 0.78692, 0.171778, 200.0, 0.990413, 0.816451, 0.00729848]
    native_fieldsdark_matter_densityLUT.ShowDataHistogram = 1
    native_fieldsdark_matter_densityLUT.ColorSpace = 'Lab'
    native_fieldsdark_matter_densityLUT.NanColor = [1.0, 0.0, 0.0]
    native_fieldsdark_matter_densityLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsdark_matter_density'
    native_fieldsdark_matter_densityPWF = GetOpacityTransferFunction('native_fieldsdark_matter_density')
    native_fieldsdark_matter_densityPWF.Points = [3.0, 0.0, 0.5, 0.0, 200.0, 1.0, 0.5, 0.0]
    native_fieldsdark_matter_densityPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage2Display.Representation = 'Volume'
    resampleToImage2Display.ColorArrayName = ['POINTS', 'native_fields/dark_matter_density']
    resampleToImage2Display.LookupTable = native_fieldsdark_matter_densityLUT
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
    resampleToImage2Display.ScalarOpacityFunction = native_fieldsdark_matter_densityPWF
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

    # get color legend/bar for native_fieldsdark_matter_densityLUT in view renderView2
    native_fieldsdark_matter_densityLUTColorBar = GetScalarBar(native_fieldsdark_matter_densityLUT, renderView2)
    native_fieldsdark_matter_densityLUTColorBar.Title = 'native_fields/dark_matter_density'
    native_fieldsdark_matter_densityLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldsdark_matter_densityLUTColorBar.Visibility = 1

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