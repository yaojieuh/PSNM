
try: paraview.simple
except: from paraview.simple import *

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    if datadescription.GetForceOutput() == True:
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    timestep = datadescription.GetTimeStep()

    input_name = 'input'
    if (timestep % 1 == 0) :
        datadescription.GetInputDescriptionByName(input_name).AllFieldsOn()
        datadescription.GetInputDescriptionByName(input_name).GenerateMeshOn()
    else:
        datadescription.GetInputDescriptionByName(input_name).AllFieldsOff()
        datadescription.GetInputDescriptionByName(input_name).GenerateMeshOff()


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    cp_writers = []
    cp_views = []
    timestep = datadescription.GetTimeStep()

    RenderView1 = CreateView( CreateRenderView, "image_%t.png", 1, 0, 1, cp_views )
    RenderView1.ViewSize = [800,600]
    RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
    RenderView1.KeyLightAzimuth = 10.0
    RenderView1.UseTexturedBackground = 0
    RenderView1.UseLight = 1
    RenderView1.CameraPosition = [2.3084675093833082, 2.6459281240458798, 2.2934846730499037]
    RenderView1.FillLightKFRatio = 3.0
    RenderView1.Background2 = [0.0, 0.0, 0.16500000000000001]
    RenderView1.FillLightAzimuth = -10.0
    RenderView1.LODResolution = 50.0
    RenderView1.BackgroundTexture = []
    RenderView1.InteractionMode = '3D'
    RenderView1.StencilCapable = 1
    RenderView1.LightIntensity = 1.0
    RenderView1.CameraFocalPoint = [0.54306134520134663, 0.44214965793510208, 0.52273486875626163]
    RenderView1.ImageReductionFactor = 2
    RenderView1.CameraViewAngle = 30.0
    RenderView1.CameraParallelScale = 0.86264249205090571
    RenderView1.EyeAngle = 2.0
    RenderView1.HeadLightKHRatio = 3.0
    RenderView1.StereoRender = 0
    RenderView1.KeyLightIntensity = 0.75
    RenderView1.BackLightAzimuth = 110.0
    RenderView1.OrientationAxesInteractivity = 0
    RenderView1.UseInteractiveRenderingForSceenshots = 0
    RenderView1.UseOffscreenRendering = 0
    RenderView1.Background = [0.31999694819562063, 0.34000152590218968, 0.42999923704890519]
    RenderView1.UseOffscreenRenderingForScreenshots = 0
    RenderView1.NonInteractiveRenderDelay = 2
    RenderView1.CenterOfRotation = [0.498046875, 0.498046875, 0.498046875]
    RenderView1.CameraParallelProjection = 0
    RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
    RenderView1.HeadLightWarmth = 0.40000000000000002
    RenderView1.MaximumNumberOfPeels = 4
    RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
    RenderView1.StereoType = 'Red-Blue'
    RenderView1.DepthPeeling = 1
    RenderView1.BackLightKBRatio = 3.5
    RenderView1.StereoCapableWindow = 1
    RenderView1.CameraViewUp = [-0.46266759705642396, 0.75019611734053027, -0.47238171023149073]
    RenderView1.LightType = 'HeadLight'
    RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.CacheKey = 0.0
    RenderView1.UseCache = 0
    RenderView1.KeyLightElevation = 50.0
    RenderView1.CenterAxesVisibility = 0
    RenderView1.MaintainLuminance = 1
    RenderView1.StillRenderImageReductionFactor = 1
    RenderView1.BackLightWarmth = 0.5
    RenderView1.FillLightElevation = -75.0
    RenderView1.MultiSamples = 0
    RenderView1.FillLightWarmth = 0.40000000000000002
    RenderView1.AlphaBitPlanes = 1
    RenderView1.LightSwitch = 1
    RenderView1.OrientationAxesVisibility = 1
    RenderView1.CameraClippingRange = [1.5928136999702023, 5.5315633152353296]
    RenderView1.BackLightElevation = 0.0
    RenderView1.ViewTime = 0.0
    RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    RenderView1.LODThreshold = 18.199999999999999
    RenderView1.CollectGeometryThreshold = 100.0
    RenderView1.UseGradientBackground = 0
    RenderView1.KeyLightWarmth = 0.59999999999999998
    RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    test2_ = CreateProducer( datadescription, "input" )
    
    Calculator1 = Calculator( guiName="Calculator1", Function='mag(realtempx*iHat+realtempy*jHat+realtempz*kHat)', ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, AttributeMode='point_data', CoordinateResults=0 )
    
    a1_realtempx_PiecewiseFunction = CreatePiecewiseFunction( Points=[-1.5000000000002331, 0.81166660785675049, 0.5, 0.74999964237213135, 0.039327789098029697, 0.0, 0.5, 0.62068736553192139, 1.5000000000002318, 0.26166665554046631, 0.5, 0.0] )
    
    a1_Result_PiecewiseFunction = CreatePiecewiseFunction( Points=[2.3316233587421763e-15, 0.82999998331069946, 0.5, 0.88019716739654541, 0.86602540378450699, 0.0, 0.5, 0.33582091331481934, 1.7320507764816284, 0.21487997472286224, 0.5, 0.0] )
    
    a1_realtempx_PVLookupTable = GetLookupTableForArray( "realtempx", 1, Discretize=1, RGBPoints=[-1.5000000000002331, 0.0, 1.0, 1.0, -0.15000000000002389, 0.0078431372549019607, 0.0, 0.94509803921568625, -4.4408920985006262e-16, 0.17254901960784313, 0.0, 0.41568627450980394, 0.150000000000023, 0.9137254901960784, 0.0, 0.043137254901960784, 1.5000000000002318, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, ColorSpace='RGB', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_Result_PVLookupTable = GetLookupTableForArray( "Result", 1, Discretize=1, RGBPoints=[2.3316233520252134e-15, 0.0, 1.0, 1.0, 0.77942286340605649, 0.011764705882352941, 0.0, 0.93333333333333335, 0.86602540378450699, 0.20000000000000001, 0.0, 0.40000000000000002, 0.95262794416295749, 0.89803921568627454, 0.0, 0.050980392156862744, 1.7320508075690118, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, ColorSpace='RGB', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    SetActiveSource(test2_)
    DataRepresentation1 = Show()
    DataRepresentation1.CubeAxesZAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation1.SelectionPointFieldDataArrayName = 'realtempx'
    DataRepresentation1.SuppressLOD = 0
    DataRepresentation1.CubeAxesXGridLines = 0
    DataRepresentation1.CubeAxesYAxisTickVisibility = 1
    DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Position = [0.0, 0.0, 0.0]
    DataRepresentation1.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation1.SelectionOpacity = 1.0
    DataRepresentation1.SelectionPointLabelShadow = 0
    DataRepresentation1.CubeAxesYGridLines = 0
    DataRepresentation1.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation1.OrientationMode = 'Direction'
    DataRepresentation1.Source.TipResolution = 6
    DataRepresentation1.ScaleMode = 'No Data Scaling Off'
    DataRepresentation1.Diffuse = 1.0
    DataRepresentation1.SelectionUseOutline = 0
    DataRepresentation1.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation1.Specular = 0.10000000000000001
    DataRepresentation1.SelectionVisibility = 1
    DataRepresentation1.InterpolateScalarsBeforeMapping = 1
    DataRepresentation1.CubeAxesZAxisTickVisibility = 1
    DataRepresentation1.Origin = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesVisibility = 0
    DataRepresentation1.Scale = [1.0, 1.0, 1.0]
    DataRepresentation1.SelectionCellLabelJustification = 'Left'
    DataRepresentation1.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Shade = 0
    DataRepresentation1.SelectionCellLabelOpacity = 1.0
    DataRepresentation1.CubeAxesInertia = 1
    DataRepresentation1.Source = "Arrow"
    DataRepresentation1.Source.Invert = 0
    DataRepresentation1.Masking = 0
    DataRepresentation1.Opacity = 1.0
    DataRepresentation1.LineWidth = 1.0
    DataRepresentation1.MeshVisibility = 0
    DataRepresentation1.Visibility = 0
    DataRepresentation1.SelectionCellLabelFontSize = 18
    DataRepresentation1.CubeAxesCornerOffset = 0.0
    DataRepresentation1.SelectionPointLabelJustification = 'Left'
    DataRepresentation1.SelectionPointLabelVisibility = 0
    DataRepresentation1.SelectOrientationVectors = ''
    DataRepresentation1.CubeAxesTickLocation = 'Inside'
    DataRepresentation1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.CubeAxesYAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation1.Source.ShaftResolution = 6
    DataRepresentation1.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation1.SelectScaleArray = ''
    DataRepresentation1.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation1.ColorAttributeType = 'POINT_DATA'
    DataRepresentation1.SpecularPower = 100.0
    DataRepresentation1.Texture = []
    DataRepresentation1.SelectionCellLabelShadow = 0
    DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.MapScalars = 1
    DataRepresentation1.PointSize = 2.0
    DataRepresentation1.Source.TipLength = 0.34999999999999998
    DataRepresentation1.SelectionCellLabelFormat = ''
    DataRepresentation1.Scaling = 0
    DataRepresentation1.StaticMode = 0
    DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation1.SliceMode = 'XY Plane'
    DataRepresentation1.Source.TipRadius = 0.10000000000000001
    DataRepresentation1.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation1.CubeAxesXAxisTickVisibility = 1
    DataRepresentation1.SelectionCellLabelVisibility = 0
    DataRepresentation1.NonlinearSubdivisionLevel = 1
    DataRepresentation1.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation1.Representation = 'Outline'
    DataRepresentation1.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation1.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation1.CubeAxesXTitle = 'X-Axis'
    DataRepresentation1.ScalarOpacityUnitDistance = 0.0067658234670659291
    DataRepresentation1.BackfaceOpacity = 1.0
    DataRepresentation1.SelectionPointLabelFontSize = 18
    DataRepresentation1.SelectionCellFieldDataArrayName = 'vtkOriginalCellIds'
    DataRepresentation1.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation1.Ambient = 0.0
    DataRepresentation1.VolumeRenderingMode = 'Smart'
    DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation1.ScaleFactor = 0.099609375
    DataRepresentation1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Slice = 0
    DataRepresentation1.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation1.ScalarOpacityFunction = []
    DataRepresentation1.SelectMaskArray = ''
    DataRepresentation1.SelectionLineWidth = 2.0
    DataRepresentation1.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation1.CubeAxesXAxisVisibility = 1
    DataRepresentation1.Interpolation = 'Gouraud'
    DataRepresentation1.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation1.SelectionCellLabelItalic = 0
    DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation1.InterpolationType = 'Linear'
    DataRepresentation1.CubeAxesZGridLines = 0
    DataRepresentation1.SelectionPointLabelFormat = ''
    DataRepresentation1.SelectionPointLabelOpacity = 1.0
    DataRepresentation1.Pickable = 1
    DataRepresentation1.CustomBoundsActive = [0, 0, 0]
    DataRepresentation1.SelectionRepresentation = 'Wireframe'
    DataRepresentation1.SelectionPointLabelBold = 0
    DataRepresentation1.ColorArrayName = ''
    DataRepresentation1.SelectionPointLabelItalic = 0
    DataRepresentation1.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation1.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation1.LookupTable = []
    DataRepresentation1.SelectionPointSize = 5.0
    DataRepresentation1.SelectionCellLabelBold = 0
    DataRepresentation1.Orient = 0
    
    SetActiveSource(Calculator1)
    DataRepresentation2 = Show()
    DataRepresentation2.CubeAxesZAxisVisibility = 1
    DataRepresentation2.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation2.SelectionPointFieldDataArrayName = 'Result'
    DataRepresentation2.SuppressLOD = 0
    DataRepresentation2.CubeAxesXGridLines = 0
    DataRepresentation2.CubeAxesYAxisTickVisibility = 1
    DataRepresentation2.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Position = [0.0, 0.0, 0.0]
    DataRepresentation2.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation2.SelectionOpacity = 1.0
    DataRepresentation2.SelectionPointLabelShadow = 0
    DataRepresentation2.CubeAxesYGridLines = 0
    DataRepresentation2.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation2.OrientationMode = 'Direction'
    DataRepresentation2.Source.TipResolution = 6
    DataRepresentation2.ScaleMode = 'No Data Scaling Off'
    DataRepresentation2.Diffuse = 1.0
    DataRepresentation2.SelectionUseOutline = 0
    DataRepresentation2.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation2.Specular = 0.10000000000000001
    DataRepresentation2.SelectionVisibility = 1
    DataRepresentation2.InterpolateScalarsBeforeMapping = 1
    DataRepresentation2.CubeAxesZAxisTickVisibility = 1
    DataRepresentation2.Origin = [0.0, 0.0, 0.0]
    DataRepresentation2.CubeAxesVisibility = 0
    DataRepresentation2.Scale = [1.0, 1.0, 1.0]
    DataRepresentation2.SelectionCellLabelJustification = 'Left'
    DataRepresentation2.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Shade = 0
    DataRepresentation2.SelectionCellLabelOpacity = 1.0
    DataRepresentation2.CubeAxesInertia = 1
    DataRepresentation2.Source = "Arrow"
    DataRepresentation2.Source.Invert = 0
    DataRepresentation2.Masking = 0
    DataRepresentation2.Opacity = 1.0
    DataRepresentation2.LineWidth = 1.0
    DataRepresentation2.MeshVisibility = 0
    DataRepresentation2.Visibility = 1
    DataRepresentation2.SelectionCellLabelFontSize = 18
    DataRepresentation2.CubeAxesCornerOffset = 0.0
    DataRepresentation2.SelectionPointLabelJustification = 'Left'
    DataRepresentation2.SelectionPointLabelVisibility = 0
    DataRepresentation2.SelectOrientationVectors = ''
    DataRepresentation2.CubeAxesTickLocation = 'Inside'
    DataRepresentation2.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation2.CubeAxesYAxisVisibility = 1
    DataRepresentation2.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation2.Source.ShaftResolution = 6
    DataRepresentation2.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation2.SelectScaleArray = ''
    DataRepresentation2.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation2.ColorAttributeType = 'POINT_DATA'
    DataRepresentation2.SpecularPower = 100.0
    DataRepresentation2.Texture = []
    DataRepresentation2.SelectionCellLabelShadow = 0
    DataRepresentation2.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation2.MapScalars = 1
    DataRepresentation2.PointSize = 2.0
    DataRepresentation2.Source.TipLength = 0.34999999999999998
    DataRepresentation2.SelectionCellLabelFormat = ''
    DataRepresentation2.Scaling = 0
    DataRepresentation2.StaticMode = 0
    DataRepresentation2.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation2.SliceMode = 'XY Plane'
    DataRepresentation2.Source.TipRadius = 0.10000000000000001
    DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation2.CubeAxesXAxisTickVisibility = 1
    DataRepresentation2.SelectionCellLabelVisibility = 0
    DataRepresentation2.NonlinearSubdivisionLevel = 1
    DataRepresentation2.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation2.Representation = 'Volume'
    DataRepresentation2.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation2.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation2.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation2.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation2.CubeAxesXTitle = 'X-Axis'
    DataRepresentation2.ScalarOpacityUnitDistance = 0.0067658234670659291
    DataRepresentation2.BackfaceOpacity = 1.0
    DataRepresentation2.SelectionPointLabelFontSize = 18
    DataRepresentation2.SelectionCellFieldDataArrayName = 'vtkOriginalCellIds'
    DataRepresentation2.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation2.Ambient = 0.0
    DataRepresentation2.VolumeRenderingMode = 'Smart'
    DataRepresentation2.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation2.ScaleFactor = 0.099609375
    DataRepresentation2.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Slice = 0
    DataRepresentation2.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation2.ScalarOpacityFunction = a1_Result_PiecewiseFunction
    DataRepresentation2.SelectMaskArray = ''
    DataRepresentation2.SelectionLineWidth = 2.0
    DataRepresentation2.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation2.CubeAxesXAxisVisibility = 1
    DataRepresentation2.Interpolation = 'Gouraud'
    DataRepresentation2.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation2.SelectionCellLabelItalic = 0
    DataRepresentation2.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation2.InterpolationType = 'Linear'
    DataRepresentation2.CubeAxesZGridLines = 0
    DataRepresentation2.SelectionPointLabelFormat = ''
    DataRepresentation2.SelectionPointLabelOpacity = 1.0
    DataRepresentation2.Pickable = 1
    DataRepresentation2.CustomBoundsActive = [0, 0, 0]
    DataRepresentation2.SelectionRepresentation = 'Wireframe'
    DataRepresentation2.SelectionPointLabelBold = 0
    DataRepresentation2.ColorArrayName = 'Result'
    DataRepresentation2.SelectionPointLabelItalic = 0
    DataRepresentation2.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation2.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation2.LookupTable = a1_Result_PVLookupTable
    DataRepresentation2.SelectionPointSize = 5.0
    DataRepresentation2.SelectionCellLabelBold = 0
    DataRepresentation2.Orient = 0
    

    for writer in cp_writers:
        if timestep % writer.cpFrequency == 0 or datadescription.GetForceOutput() == True:
            writer.FileName = writer.cpFileName.replace("%t", str(timestep))
            writer.UpdatePipeline()

    if False : # rescale data range
        import math
        for view in cp_views:
            if timestep % view.cpFrequency == 0 or datadescription.GetForceOutput() == True:
                reps = view.Representations
                for rep in reps:
                    if hasattr(rep, 'Visibility') and rep.Visibility == 1 and hasattr(rep, 'MapScalars') and rep.MapScalars != '':
                        input = rep.Input
                        input.UpdatePipeline() #make sure range is up-to-date
                        lut = rep.LookupTable
                        if lut == None:
                            continue
                        if rep.ColorAttributeType == 'POINT_DATA':
                            datainformation = input.GetPointDataInformation()
                        elif rep.ColorAttributeType == 'CELL_DATA':
                            datainformation = input.GetCellDataInformation()
                        else:
                            print 'something strange with color attribute type', rep.ColorAttributeType

                        if lut.VectorMode != 'Magnitude' or                            datainformation.GetArray(rep.ColorArrayName).GetNumberOfComponents() == 1:
                            datarange = datainformation.GetArray(rep.ColorArrayName).GetRange(lut.VectorComponent)
                        else:
                            datarange = [0,0]
                            for i in range(datainformation.GetArray(rep.ColorArrayName).GetNumberOfComponents()):
                                for j in range(2):
                                    datarange[j] += datainformation.GetArray(rep.ColorArrayName).GetRange(i)[j]*datainformation.GetArray(rep.ColorArrayName).GetRange(i)[j]
                            datarange[0] = math.sqrt(datarange[0])
                            datarange[1] = math.sqrt(datarange[1])

                        rgbpoints = lut.RGBPoints.GetData()
                        numpts = len(rgbpoints)/4
                        minvalue = min(datarange[0], rgbpoints[0])
                        maxvalue = max(datarange[1], rgbpoints[(numpts-1)*4])
                        if minvalue != rgbpoints[0] or maxvalue != rgbpoints[(numpts-1)*4]:
                            # rescale all of the points
                            oldrange = rgbpoints[(numpts-1)*4] - rgbpoints[0]
                            newrange = maxvalue - minvalue
                            newrgbpoints = list(rgbpoints)
                            for v in range(numpts):
                                newrgbpoints[v*4] = minvalue+(rgbpoints[v*4] - rgbpoints[0])*newrange/oldrange

                            lut.RGBPoints.SetData(newrgbpoints)

    for view in cp_views:
        if timestep % view.cpFrequency == 0 or datadescription.GetForceOutput() == True:
            fname = view.cpFileName
            fname = fname.replace("%t", str(timestep))
            if view.cpFitToScreen != 0:
                if view.IsA("vtkSMRenderViewProxy") == True:
                    view.ResetCamera()
                elif view.IsA("vtkSMContextViewProxy") == True:
                    view.ResetDisplay()
                else:
                    print ' do not know what to do with a ', view.GetClassName()

            WriteImage(fname, view, Magnification=view.cpMagnification)


    # explicitly delete the proxies -- we do it this way to avoid problems with prototypes
    tobedeleted = GetNextProxyToDelete()
    while tobedeleted != None:
        Delete(tobedeleted)
        tobedeleted = GetNextProxyToDelete()

def GetNextProxyToDelete():
    proxyiterator = servermanager.ProxyIterator()
    for proxy in proxyiterator:
        group = proxyiterator.GetGroup()
        if group.find("prototypes") != -1:
            continue
        if group != 'timekeeper' and group.find("pq_helper_proxies") == -1 :
            return proxy
    return None

def CreateProducer(datadescription, gridname):
    "Creates a producer proxy for the grid"
    if not datadescription.GetInputDescriptionByName(gridname):
        raise RuntimeError, "Simulation input name '%s' does not exist" % gridname
    grid = datadescription.GetInputDescriptionByName(gridname).GetGrid()
    producer = PVTrivialProducer()
    producer.GetClientSideObject().SetOutput(grid)
    if grid.IsA("vtkImageData") == True or grid.IsA("vtkStructuredGrid") == True or grid.IsA("vtkRectilinearGrid") == True:
        extent = datadescription.GetInputDescriptionByName(gridname).GetWholeExtent()
        producer.WholeExtent= [ extent[0], extent[1], extent[2], extent[3], extent[4], extent[5] ]

    producer.UpdatePipeline()
    return producer


def CreateWriter(proxy_ctor, filename, freq, cp_writers):
    writer = proxy_ctor()
    writer.FileName = filename
    writer.add_attribute("cpFrequency", freq)
    writer.add_attribute("cpFileName", filename)
    cp_writers.append(writer)
    return writer

def CreateView(proxy_ctor, filename, freq, fittoscreen, magnification, cp_views):
    view = proxy_ctor()
    view.add_attribute("cpFileName", filename)
    view.add_attribute("cpFrequency", freq)
    view.add_attribute("cpFileName", filename)
    view.add_attribute("cpFitToScreen", fittoscreen)
    view.add_attribute("cpMagnification", magnification)
    cp_views.append(view)
    return view
