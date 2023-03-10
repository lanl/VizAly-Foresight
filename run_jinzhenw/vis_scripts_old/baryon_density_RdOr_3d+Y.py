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
    nVB_C009_l10n512_S12345T692_z54h5.CellArrays = ['native_fields/baryon_density']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=nVB_C009_l10n512_S12345T692_z54h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityLUT = GetColorTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityLUT.AutomaticRescaleRangeMode = 'Never'
    native_fieldsbaryon_densityLUT.RGBPoints = [0.06377197802066803, 0.498039, 0.0, 0.0, 31.43227060525976, 0.6004, 0.0, 0.0, 62.80076923249885, 0.702514, 0.000738, 0.000477, 94.16951782785196, 0.773379, 0.095225, 0.061499, 125.53801645509105, 0.843875, 0.189865, 0.12283, 156.90651508233015, 0.891119, 0.294195, 0.203537, 188.27501370956924, 0.937855, 0.397924, 0.283137, 219.64351233680833, 0.963445, 0.476663, 0.316601, 251.0121534458724, 0.988297, 0.555771, 0.351665, 282.38075955940053, 0.990265, 0.646321, 0.436309, 313.7492581866396, 0.992157, 0.735256, 0.519646, 345.1177568138787, 0.992157, 0.784468, 0.570827, 376.48625544111786, 0.992249, 0.833218, 0.623483, 407.8550040364709, 0.994218, 0.872587, 0.706159, 439.22350266371, 0.996186, 0.911419, 0.788189, 470.59200129094916, 0.998155, 0.940946, 0.859054, 500.0, 1.0, 0.968627, 0.92549]
    native_fieldsbaryon_densityLUT.ShowDataHistogram = 1
    native_fieldsbaryon_densityLUT.AutomaticDataHistogramComputation = 1
    native_fieldsbaryon_densityLUT.ColorSpace = 'Lab'
    native_fieldsbaryon_densityLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityPWF = GetOpacityTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityPWF.Points = [0.06377197802066803, 0.0, 0.5, 0.0, 14.26650619506836, 0.13368983566761017, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]
    native_fieldsbaryon_densityPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Volume'
    resampleToImage1Display.ColorArrayName = ['POINTS', 'native_fields/baryon_density']
    resampleToImage1Display.LookupTable = native_fieldsbaryon_densityLUT
    resampleToImage1Display.SelectTCoordArray = 'None'
    resampleToImage1Display.SelectNormalArray = 'None'
    resampleToImage1Display.SelectTangentArray = 'None'
    resampleToImage1Display.OSPRayScaleArray = 'native_fields/baryon_density'
    resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleToImage1Display.SelectOrientationVectors = 'None'
    resampleToImage1Display.ScaleFactor = 51.199948799999994
    resampleToImage1Display.SelectScaleArray = 'None'
    resampleToImage1Display.GlyphType = 'Arrow'
    resampleToImage1Display.GlyphTableIndexArray = 'None'
    resampleToImage1Display.GaussianRadius = 2.5599974399999996
    resampleToImage1Display.SetScaleArray = ['POINTS', 'native_fields/baryon_density']
    resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.OpacityArray = ['POINTS', 'native_fields/baryon_density']
    resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
    resampleToImage1Display.ScalarOpacityUnitDistance = 1.135438604041588
    resampleToImage1Display.ScalarOpacityFunction = native_fieldsbaryon_densityPWF
    resampleToImage1Display.OpacityArrayName = ['POINTS', 'native_fields/baryon_density']
    resampleToImage1Display.SliceFunction = 'Plane'
    resampleToImage1Display.Slice = 255

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage1Display.ScaleTransferFunction.Points = [0.06377197802066803, 0.0, 0.5, 0.0, 28918.84765625, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage1Display.OpacityTransferFunction.Points = [0.06377197802066803, 0.0, 0.5, 0.0, 28918.84765625, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    resampleToImage1Display.SliceFunction.Origin = [256.0, 256.0, 256.0]

    # setup the color legend parameters for each legend in this view

    # get color legend/bar for native_fieldsbaryon_densityLUT in view renderView1
    native_fieldsbaryon_densityLUTColorBar = GetScalarBar(native_fieldsbaryon_densityLUT, renderView1)
    native_fieldsbaryon_densityLUTColorBar.Title = 'native_fields/baryon_density'
    native_fieldsbaryon_densityLUTColorBar.ComponentTitle = ''

    # set color bar visibility
    native_fieldsbaryon_densityLUTColorBar.Visibility = 1

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
