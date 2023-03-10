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

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # Create a new 'Render View'
    renderView1 = CreateView('RenderView')
    renderView1.ViewSize = [1424, 955]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.CenterOfRotation = [256.0, 256.0, 256.0]
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraPosition = [-1457.1836768696405, 256.0, 256.0]
    renderView1.CameraFocalPoint = [256.0, 256.0, 256.0]
    renderView1.CameraViewUp = [-0.31088313855757344, 0.7211703487654267, 0.619084002557121]
    renderView1.CameraFocalDisk = 1.0
    renderView1.CameraParallelScale = 443.40500673763256
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.OSPRayMaterialLibrary = materialLibrary1

    SetActiveView(None)

    colorPalette = GetSettingsProxy('ColorPalette')
    colorPalette.Background = [0.0, 0.0, 0.0]

    # ----------------------------------------------------------------
    # setup view layouts
    # ----------------------------------------------------------------

    # create new layout object 'Layout #1'
    layout1 = CreateLayout(name='Layout #1')
    layout1.AssignView(0, renderView1)
    layout1.SetSize(1424, 955)

    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(renderView1)
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------

    # create a new 'VisItPixieReader'
    nVB_C009_l10n512_S12345T692_z42h5 = VisItPixieReader(registrationName=regName, FileName=inputFilename)
    nVB_C009_l10n512_S12345T692_z42h5.Meshes = ['mesh_512x512x512']
    nVB_C009_l10n512_S12345T692_z42h5.CellArrays = ['native_fields/baryon_density']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=nVB_C009_l10n512_S12345T692_z42h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityLUT = GetColorTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityLUT.RGBPoints = [0.04938583821058273, 0.0, 1.0, 1.0, 1106.7252197265625, 0.0, 0.9294117647058824, 1.0, 2471.625244140625, 0.0, 0.0, 1.0, 3578.30126953125, 0.0, 0.0, 0.501960784314, 4943.201171875, 1.0, 0.0, 0.0, 6418.76904296875, 1.0, 0.9333333333333333, 0.0, 71011.75, 1.0, 1.0, 0.0]
    native_fieldsbaryon_densityLUT.ShowDataHistogram = 1
    native_fieldsbaryon_densityLUT.ColorSpace = 'RGB'
    native_fieldsbaryon_densityLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsbaryon_density'
    native_fieldsbaryon_densityPWF = GetOpacityTransferFunction('native_fieldsbaryon_density')
    native_fieldsbaryon_densityPWF.Points = [0.04938583821058273, 0.0, 0.5, 0.0, 548.7048950195312, 0.9518716931343079, 0.5, 0.0, 71011.75, 1.0, 0.5, 0.0]
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
    resampleToImage1Display.ScalarOpacityUnitDistance = 0.3735438604041588
    resampleToImage1Display.ScalarOpacityFunction = native_fieldsbaryon_densityPWF
    resampleToImage1Display.OpacityArrayName = ['POINTS', 'native_fields/baryon_density']
    resampleToImage1Display.SliceFunction = 'Plane'
    resampleToImage1Display.Slice = 49

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage1Display.ScaleTransferFunction.Points = [0.05396782234311104, 0.0, 0.5, 0.0, 45861.0703125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage1Display.OpacityTransferFunction.Points = [0.05396782234311104, 0.0, 0.5, 0.0, 45861.0703125, 1.0, 0.5, 0.0]

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

    SaveScreenshot(outputName, renderView1, ImageResolution=[716, 680], CompressionLevel='0')

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
