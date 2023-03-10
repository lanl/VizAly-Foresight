# state file generated using paraview version 5.10.0
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
    renderView1.CameraPosition = [857.5290395944288, 1230.2583713062552, -576.8449960429356]
    renderView1.CameraFocalPoint = [256.00000000000006, 255.99999999999963, 255.99999999999991]
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
    #nVB_C009_l10n512_S12345T692_z42h5 = VisItPixieReader(registrationName='NVB_C009_l10n512_S12345T692_z42.h5', FileName='/home/pascal/data/NVB_C009_l10n512_S12345T692_z42.h5')
    nVB_C009_l10n512_S12345T692_z42h5 = VisItPixieReader(registrationName=regName, FileName=inputFilename)
    nVB_C009_l10n512_S12345T692_z42h5.Meshes = ['mesh_512x512x512']
    nVB_C009_l10n512_S12345T692_z42h5.CellArrays = ['native_fields/dark_matter_density']

    # create a new 'Resample To Image'
    resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=nVB_C009_l10n512_S12345T692_z42h5)
    resampleToImage1.SamplingDimensions = [512, 512, 512]
    resampleToImage1.SamplingBounds = [0.0, 512.0, 0.0, 512.0, 0.0, 512.0]

    # ----------------------------------------------------------------
    # setup the visualization in view 'renderView1'
    # ----------------------------------------------------------------

    # show data from resampleToImage1
    resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

    # get color transfer function/color map for 'native_fieldsdark_matter_density'
    native_fieldsdark_matter_densityLUT = GetColorTransferFunction('native_fieldsdark_matter_density')
    native_fieldsdark_matter_densityLUT.EnableOpacityMapping = 1
    native_fieldsdark_matter_densityLUT.RGBPoints = [0.0, 0.0, 1.0, 1.0, 3513.2209716796874, 0.0, 0.0, 1.0, 3903.578857421875, 0.0, 0.0, 0.501960784314, 4293.936743164063, 1.0, 0.0, 0.0, 7807.15771484375, 1.0, 1.0, 0.0]
    native_fieldsdark_matter_densityLUT.ShowDataHistogram = 1
    native_fieldsdark_matter_densityLUT.ColorSpace = 'RGB'
    native_fieldsdark_matter_densityLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'native_fieldsdark_matter_density'
    native_fieldsdark_matter_densityPWF = GetOpacityTransferFunction('native_fieldsdark_matter_density')
    native_fieldsdark_matter_densityPWF.Points = [0.0, 0.0, 0.5, 0.0, 70.5466079711914, 1.0, 0.5, 0.0, 7807.15771484375, 1.0, 0.5, 0.0]
    native_fieldsdark_matter_densityPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    resampleToImage1Display.Representation = 'Volume'
    resampleToImage1Display.ColorArrayName = ['POINTS', 'native_fields/dark_matter_density']
    resampleToImage1Display.LookupTable = native_fieldsdark_matter_densityLUT
    resampleToImage1Display.SelectTCoordArray = 'None'
    resampleToImage1Display.SelectNormalArray = 'None'
    resampleToImage1Display.SelectTangentArray = 'None'
    resampleToImage1Display.OSPRayScaleArray = 'native_fields/dark_matter_density'
    resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    resampleToImage1Display.SelectOrientationVectors = 'None'
    resampleToImage1Display.ScaleFactor = 51.199948799999994
    resampleToImage1Display.SelectScaleArray = 'None'
    resampleToImage1Display.GlyphType = 'Arrow'
    resampleToImage1Display.GlyphTableIndexArray = 'None'
    resampleToImage1Display.GaussianRadius = 2.5599974399999996
    resampleToImage1Display.SetScaleArray = ['POINTS', 'native_fields/dark_matter_density']
    resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.OpacityArray = ['POINTS', 'native_fields/dark_matter_density']
    resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
    resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
    resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
    resampleToImage1Display.ScalarOpacityUnitDistance = 1.7354386040415881
    resampleToImage1Display.ScalarOpacityFunction = native_fieldsdark_matter_densityPWF
    resampleToImage1Display.OpacityArrayName = ['POINTS', 'native_fields/dark_matter_density']
    resampleToImage1Display.SliceFunction = 'Plane'
    resampleToImage1Display.Slice = 255

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    resampleToImage1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7807.15771484375, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    resampleToImage1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 7807.15771484375, 1.0, 0.5, 0.0]

    # init the 'Plane' selected for 'SliceFunction'
    resampleToImage1Display.SliceFunction.Origin = [256.0, 256.0, 256.0]

    # save screenshot
    #SaveScreenshot('/home/pascal/Desktop/test3.png', renderView1, ImageResolution=[716, 680], CompressionLevel='0')
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

# Run as:
# ../ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin/pvpython nyx_img_test.py /home/pascal/data/NVB_C009_l10n512_S12345T692_z42.h5 NVB_C009_l10n512_S12345T692_z42.h5 /home/pascal/Desktop/test4.png

# View Image
# ssh -Y darwin-fe
# display test.png


