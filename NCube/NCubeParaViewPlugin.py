# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import *
#from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase

def _str(text):
    import sys
    # fix string issue for Python 2
    if sys.version_info < (3, 0) and hasattr(text, 'encode'):
        return text.encode('utf-8')
    return str(text)

# Load shapefile or geojson
def _NCubeGeoDataFrameLoad(shapename, shapecol=None, shapeencoding=None, extent=None, dem_crs=None):
    import geopandas as gpd
    from functools import partial
    from shapely.ops import transform

    df = gpd.read_file(shapename, encoding=shapeencoding)
    # remove NULL geometries
    df = df[df.geometry.notnull()]
    if (len(df)) == 0:
        return
    if shapecol is not None:
        df = df.sort_values(shapecol).set_index(shapecol)
    #print ("shapecol",shapecol)

    # extract the geometry coordinate system
    df_crs = str(df.crs['init']) if df.crs != {} else None
    print ("df_crs",df_crs,"dem_crs",dem_crs)

    # reproject when the both coordinate systems are defined and these are different
    if df_crs and dem_crs:
        df_extent = gpd.GeoDataFrame([], crs={'init' : dem_crs}, geometry=[extent])
        extent_reproj = df_extent.to_crs({'init' : df_crs})['geometry'][0]
        # if original or reprojected raster extent is valid, use it to crop geometry
        if extent_reproj.is_valid:
            # geometry intersection to raster extent in geometry coordinate system
            df = df[df.geometry.intersects(extent_reproj)]
            df['geometry'] = df.geometry.intersection(extent_reproj)

        # reproject [cropped] geometry to original raster coordinates if needed
        return df.to_crs({'init' : dem_crs})

    # let's assume the coordinate systems are the same
    if extent is not None:
        df = df[df.geometry.intersects(extent)]
        df['geometry'] = df.geometry.intersection(extent)

    return df

def _NcubeDataFrameToVTKArrays(df):
    from vtk import vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray

    arrays = []
    # Create columns
    for colname in df.columns:
        dtype = df[colname].dtype
        #print (colname, dtype)
        if dtype in ['O','str','datetime64']:
            vtk_arr = vtkStringArray()
        elif dtype in ['int64']:
            vtk_arr = vtkIntArray()
        elif dtype in ['float64']:
            vtk_arr = vtkFloatArray()
        elif dtype in ['bool']:
            vtk_arr = vtkBitArray()
        else:
            print ('Unknown Pandas column type', dtype)
            vtk_arr = vtkStringArray()
        vtk_arr.SetNumberOfComponents(1)
        vtk_arr.SetName(colname)
        for val in df[colname]:
            # some different datatypes could be saved as strings
            if isinstance(vtk_arr, vtkStringArray):
                val = str(val)
            vtk_arr.InsertNextValue(val)
        arrays.append(vtk_arr)

    return arrays

# list of list of VtkArray's
def _NCubeGeoDataFrameRowToVTKArrays(items):
    #vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, 
    from vtk import vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry

    vtk_row = []
    for (key,value) in items.items():
        #print (key,value)
        # define attribute as array
        if isinstance(value, (BaseMultipartGeometry)):
            #print ('BaseMultipartGeometry')
            continue
        elif isinstance(value, (BaseGeometry)):
            #print ('BaseGeometry')
            continue
        elif isinstance(value, (int)):
            vtk_arr = vtkIntArray()
        elif isinstance(value, (float)):
            vtk_arr = vtkFloatArray()
        elif isinstance(value, (bool)):
            vtk_arr = vtkBitArray()
        else:
            # some different datatypes could be saved as strings
            value = _str(value)
            vtk_arr = vtkStringArray()

        vtk_arr.SetNumberOfComponents(1)
        vtk_arr.SetName(key)
        vtk_row.append((vtk_arr, value))
    return vtk_row

# process [multi]geometry
def _NCubeGeometryToPolyData(geometry, dem=None):
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    import xarray as xr
    import numpy as np

    vtk_points = vtkPoints()
    vtk_cells = vtkCellArray()
    # get part(s) of (multi)geometry
    if isinstance(geometry, (BaseMultipartGeometry)):
        geometries = [geom for geom in geometry]
    else:
        geometries = [geometry]
    print ([geom.type for geom in geometries])
    for geom in geometries:
        coords = geom.coords
        #coords = geom.exterior.coords
        xs = np.array(coords.xy[0])
        ys = np.array(coords.xy[1])
        if len(coords.xy) > 2:
            zs = np.array(coords.xy[2])
        else:
            zs = np.array([0]*len(xs))
        #print (xs)
        # rasterize geometries as lines
        if dem is not None:
#            print (dem)
            zs = dem.sel(x=xr.DataArray(xs), y=xr.DataArray(ys), method='nearest').values
        #print ("xs", xs)
        mask = np.where(~np.isnan(zs))[0]
        mask2 = np.where(np.diff(mask)!=1)[0]+1
        xs = np.split(xs[mask], mask2)
        ys = np.split(ys[mask], mask2)
        zs = np.split(zs[mask], mask2)
        for (_xs,_ys,_zs) in zip(xs,ys,zs):
            # need to have 2 point or more
            if len(_xs) <= 1:
                continue
            vtk_cells.InsertNextCell(len(_xs))
            for (x,y,z) in zip(_xs,_ys,_zs):
                pointId = vtk_points.InsertNextPoint(x, y, z)
                vtk_cells.InsertCellPoint(pointId)

    # not enougth valid points
    if vtk_points.GetNumberOfPoints() < 2:
        return

    print ("GetNumberOfPoints", vtk_points.GetNumberOfPoints())
    vtk_polyData = vtkPolyData()
    vtk_polyData.SetPoints(vtk_points)
    if geometry.type in ['Point','MultiPoint']:
        vtk_polyData.SetVerts(vtk_cells)
    else:
        vtk_polyData.SetLines(vtk_cells)

    return vtk_polyData


def _NCubeGeometryOnTopography(shapename, toponame, shapecol, shapeencoding):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    import xarray as xr
    import numpy as np

    #print ("_NCUBEGeometryOnTopography start")

    # Load shapefile
    if shapename is None:
        return

    # load DEM
    dem = dem_extent = dem_crs = None
    if toponame is not None:
        #dem = xr.open_rasterio(toponame, chunks=10000000).squeeze()
        dem = xr.open_rasterio(toponame).squeeze()
        # dask array can't be processed by this way
        dem.values[dem.values == dem.nodatavals[0]] = np.nan

        # NaN border to easy lookup
        dem.values[0,:]  = np.nan
        dem.values[-1,:] = np.nan
        dem.values[:,0]  = np.nan
        dem.values[:,-1] = np.nan

        dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())
        dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None

        #print (dem.values)

    df = _NCubeGeoDataFrameLoad(shapename, shapecol, shapeencoding, dem_extent, dem_crs)
    if df is None:
        return

    groups = df.index.unique() ;#[11454:11455]
    #print ("groups",groups)

    # TEST
#    groups = groups[:1]

    # iterate blocks
    vtk_blocks = []
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
        if hasattr(group, 'encode'):
            # select only equals
            _df = df[df.index.str.startswith(group)&df.index.str.endswith(group)].reset_index()
        else:
            _df = df[df.index == group].reset_index()
        vtk_appendPolyData = vtkAppendPolyData()
        # iterate rows with the same attributes and maybe multiple geometries
        for rowidx,row in _df.iterrows():
            #geoms = _NCubeTopographyCheckGeometries([row.geometry])
            vtk_polyData = _NCubeGeometryToPolyData(row.geometry, dem)
            if vtk_polyData is None:
                continue
            vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays(row.to_dict())
            for (vtk_arr, val) in vtk_arrays:
                for _ in range(vtk_polyData.GetNumberOfCells()):
                    vtk_arr.InsertNextValue(val)
                vtk_polyData.GetCellData().AddArray(vtk_arr)
            # compose vtkPolyData
            vtk_appendPolyData.AddInputData(vtk_polyData)
        # nothing to process
        if vtk_appendPolyData.GetNumberOfInputConnections(0) == 0:
            continue
        vtk_appendPolyData.Update()
        vtk_block = vtk_appendPolyData.GetOutput()

        vtk_blocks.append((_str(group),vtk_block))

    #print ("_NCUBEGeometryOnTopography end")

    return vtk_blocks


# TODO
# df.loc[0,'geometry']
def _NCubeTopographyToGrid(dem):
    from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, VTK_FLOAT
    from vtk.util import numpy_support as vn
    import numpy as np

    xs = dem.x.values
    ys = dem.y.values
    if dem.dims == ('y','x'):
        values = dem.values
    else:
        values = dem.values.T

    # create raster mask by geometry and for NaNs
    (yy,xx) = np.meshgrid(ys, xs)
    vtk_points = vtkPoints()
    points = np.column_stack((xx.ravel('F'),yy.ravel('F'),values.ravel('C')))
    _points = vn.numpy_to_vtk(points, deep=True)
    vtk_points.SetData(_points)
#    for (_x,_y,_z,_m) in zip(xx.ravel('F'),yy.ravel('F'),values.ravel('C'),mask.ravel('C')):
#        vtk_points.InsertNextPoint(_x,_y, _z if _m else np.nan)

    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(len(xs), len(ys), 1)
    sgrid.SetPoints(vtk_points)

    array = vn.numpy_to_vtk(values.ravel(), deep=True, array_type=VTK_FLOAT)
    array.SetName("z")
    sgrid.GetPointData().AddArray(array)

    thresh = vtkThreshold()
    thresh.SetInputData(sgrid)
    thresh.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, "z")
    thresh.ThresholdBetween(-1e30, 1e30)
    thresh.Update()

#    return sgrid
    return thresh.GetOutput()

def _NCubeTopography(shapename, toponame, shapecol, shapeencoding):
    from vtk import vtkAppendFilter
    #vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    import xarray as xr
    import numpy as np
    import rasterio
    import rasterio.mask

    print ("_NCubeTopography start")

    # Load shapefile
    if toponame is None:
        return

    # process the full topography raster
    if shapename is None:
        dem = xr.open_rasterio(toponame).squeeze()
        dem.values[dem.values == dem.nodatavals[0]] = np.nan
        vtk_ugrid = _NCubeTopographyToGrid(dem)
        print ("vtk_ugrid",vtk_ugrid)
        return [(_str('None'),vtk_ugrid)]

    # process shapefile
    dem = rasterio.open(toponame)
    print (dem.crs, dem.bounds)
    dem_crs = dem.crs.to_string() if dem.crs is not None else None
    dem_extent = box(dem.bounds.left,dem.bounds.bottom,dem.bounds.right,dem.bounds.top)
    df = _NCubeGeoDataFrameLoad(shapename, shapecol, shapeencoding, dem_extent, dem_crs)
    if df is None:
        return

    groups = df.index.unique()
    #print ("groups",groups)

    vtk_blocks = []
    # iterate blocks
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
        if hasattr(group, 'encode'):
            geoms = df[df.index.str.match(group)].geometry
        else:
            geoms = df[df.index == group].geometry

        vtk_append = vtkAppendFilter()
        # rasterize geomemetry
        try:
            out_image, out_transform = rasterio.mask.mask(dem, geoms, crop=True, filled=True)
        except:
            # geometry outside of the raster, etc.
            print (group, "exception")
            continue
        image = out_image.squeeze()
        print (image.shape)
        (xs, _) = rasterio.transform.xy(out_transform, image.shape[1]*[0], range(image.shape[1]), offset='center')
        (_, ys) = rasterio.transform.xy(out_transform, range(image.shape[0]), image.shape[0]*[0], offset='center')
        da = xr.DataArray(image, coords={'x': xs, 'y': ys}, dims=('y', 'x'))
        # process rasterized geometry
        vtk_ugrid = _NCubeTopographyToGrid(da)
        if vtk_ugrid is None:
            continue
        if shapecol:
            vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays({shapecol:group})
            for (vtk_arr, val) in vtk_arrays:
                for _ in range(vtk_ugrid.GetNumberOfCells()):
                    vtk_arr.InsertNextValue(val)
                vtk_ugrid.GetCellData().AddArray(vtk_arr)
        # compose
        vtk_append.AddInputData(vtk_ugrid)
        # nothing to process
        if vtk_append.GetNumberOfInputConnections(0) == 0:
            continue
        vtk_append.Update()
        vtk_block = vtk_append.GetOutput()
        vtk_blocks.append((_str(group),vtk_block))

    print ("_NCubeTopography end")

    return vtk_blocks


def _NCubeDataSetToGeoDataFrame(vtk_data):
        import pandas as pd
        import geopandas as gpd
        from shapely.geometry import Point, box
        from vtk.util import numpy_support
        # save all the points or only the boundaries
        if vtk_data.GetNumberOfPoints() == 0 or vtk_data.GetNumberOfPoints() > 1e6:
            # geometry is empty or it's too large (>1M points) - save only boundaries
            (minx, maxx, miny, maxy, minz, maxz) = vtk_data.GetBounds()
            geom = box(minx, miny, maxx, maxy)
            gdf = gpd.GeoDataFrame([],geometry=[geom])
        else:
            # save all points with attributes
            coords = vtk_data.GetPoints().GetData()
            ncoords = numpy_support.vtk_to_numpy(coords)
            df = pd.DataFrame(data=ncoords,columns=('x','y','z'))

            for idx in range(vtk_data.GetPointData().GetNumberOfArrays()):
                col = vtk_data.GetPointData().GetArrayName(idx)
                values = vtk_data.GetPointData().GetArray(idx)
                nvalues = numpy_support.vtk_to_numpy(values)
                df[col] = nvalues

            geom = gpd.GeoSeries(map(Point, zip(df.x, df.y, df.z)))
            df.drop(['x', 'y'], axis=1, inplace=True)
            gdf = gpd.GeoDataFrame(df, geometry=geom)

        print (gdf.head())
        return gdf


#------------------------------------------------------------------------------
# N-Cube Shapefile On Topography Block Source
#------------------------------------------------------------------------------

@smproxy.source(name="NCubeGeometryOnTopographyBlockSource",
       label="N-Cube Shapefile On Topography Block Source")
class NCubeGeometryOnTopographyBlockSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                nInputPorts=0,
                nOutputPorts=1,
                outputType='vtkMultiBlockDataSet')
        self._shapename = None
        self._shapeencoding = None
        self._shapecol = None
        self._toponame = None


    def RequestData(self, request, inInfo, outInfo):
        from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
        import time

        if self._shapename is None:
            return 1

        t0 = time.time()
        vtk_blocks = _NCubeGeometryOnTopography(self._shapename, self._toponame, self._shapecol, self._shapeencoding)
        if vtk_blocks is None or vtk_blocks == []:
            t1 = time.time()
            print ("t1-t0", t1-t0)
            return 1
        print ("vtk_blocks", len(vtk_blocks))
        mb = vtkMultiBlockDataSet.GetData(outInfo, 0)
        mb.SetNumberOfBlocks(len(vtk_blocks))
        rowidx = 0
        for (label, polyData) in vtk_blocks:
            #print (rowidx, label)
            mb.SetBlock( rowidx, polyData )
            mb.GetMetaData( rowidx ).Set( vtkCompositeDataSet.NAME(), label)
            rowidx += 1
        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1

    @smproperty.stringvector(name="Shapefile Name")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["shp", "geojson"], file_description="ESRI Shapefile, GeoJSON")
    def SetShapeFileName(self, name):
        """Specify filename for the shapefile to read."""
        print ("SetShapeFileName", name)
        name = name if name != 'None' else None
        if self._shapename != name:
            self._shapename = name
            self._shapecol = None
            self.Modified()

    @smproperty.stringvector(name="Topography File Name (optional)")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["tif", "TIF", "nc"], file_description="GeoTIFF, NetCDF")
    def SetTopographyFileName(self, name):
        """Specify filename for the topography file to read."""
        print ("SetTopographyFileName", name)
        name = name if name != 'None' else None
        if self._toponame != name:
            self._toponame = name
            self.Modified()

    @smproperty.stringvector(name="ShapeLabels", information_only="1")
    def ShapeLabels(self):
        if self._shapename is None:
            return []
        # Load shapefile
        import geopandas as gpd
        df = gpd.read_file(self._shapename, encoding=self._shapeencoding)
        cols = sorted(df.columns.values)
        del df
        if 'geometry' in cols:
            cols.remove('geometry')
        #print(cols)
        return ['None'] + list(map(str,cols))

    @smproperty.stringvector(name="Group by Field (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="ShapeLabels" function="GetShapeLabels"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetShapeLabel(self, label):
        label = label if label != 'None' else None
        self._shapecol = label
        print("SetShapeLabel ", label)
        self.Modified()




#------------------------------------------------------------------------------
# N-Cube Topography Block Source
#------------------------------------------------------------------------------

@smproxy.source(name="NCubeTopographyBlockSource",
       label="N-Cube Topography Block Source")
class NCubeTopographyBlockSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                nInputPorts=0,
                nOutputPorts=1,
                outputType='vtkMultiBlockDataSet')
        self._shapename = None
        self._shapeencoding = None
        self._shapecol = None
        self._toponame = None


    def RequestData(self, request, inInfo, outInfo):
        from vtk import vtkCompositeDataSet, vtkMultiBlockDataSet
        import time

        if self._toponame is None:
            return 1

        t0 = time.time()
        vtk_blocks = _NCubeTopography(self._shapename, self._toponame, self._shapecol, self._shapeencoding)
        if vtk_blocks is None or vtk_blocks == []:
            return 1
        print ("vtk_blocks", len(vtk_blocks))
        mb = vtkMultiBlockDataSet.GetData(outInfo, 0)
        mb.SetNumberOfBlocks(len(vtk_blocks))
        rowidx = 0
        for (label, polyData) in vtk_blocks:
            #print (rowidx, label)
            mb.SetBlock( rowidx, polyData )
            mb.GetMetaData( rowidx ).Set( vtkCompositeDataSet.NAME(), label)
            rowidx += 1
        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1

    @smproperty.stringvector(name="Shapefile Name (optional)")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["shp", "geojson"], file_description="ESRI Shapefile, GeoJSON")
    def SetShapeFileName(self, name):
        """Specify filename for the shapefile to read."""
        print ("SetShapeFileName", name)
        name = name if name != 'None' else None
        if self._shapename != name:
            self._shapename = name
            self._shapecol = None
            self.Modified()

    @smproperty.stringvector(name="Topography File Name")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["tif", "TIF", "nc"], file_description="GeoTIFF, NetCDF")
    def SetTopographyFileName(self, name):
        """Specify filename for the topography file to read."""
        print ("SetTopographyFileName", name)
        name = name if name != 'None' else None
        if self._toponame != name:
            self._toponame = name
            self.Modified()

    @smproperty.stringvector(name="ShapeLabels", information_only="1")
    def ShapeLabels(self):
        if self._shapename is None:
            return []
        # Load shapefile
        import geopandas as gpd
        df = gpd.read_file(self._shapename, encoding=self._shapeencoding)
        cols = sorted(df.columns.values)
        del df
        if 'geometry' in cols:
            cols.remove('geometry')
        #print(cols)
        return ['None'] + list(map(str,cols))

    @smproperty.stringvector(name="Group by Field (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="ShapeLabels" function="GetShapeLabels"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetShapeLabel(self, label):
        label = label if label != 'None' else None
        self._shapecol = label
        print("SetShapeLabel ", label)
        self.Modified()



#------------------------------------------------------------------------------
# N-Cube DataSet to Shapefile Writer
# N-Cube DataSet to GeoPackage Writer
#------------------------------------------------------------------------------
@smproxy.writer(extensions="shp", file_description="ESRI Shapefile", support_reload=False)
#@smproxy.writer(extensions="gpkg", file_description="GeoPackage", support_reload=False)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class NCubeShapefileWriter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=0, inputType='vtkDataSet')
        self._filename = None

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, fname):
        """Specify filename for the file to write."""
        if self._filename != fname:
            self._filename = fname
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkDataSet

        vtk_data = vtkDataSet.GetData(inInfoVec[0], 0)
        gdf = _NCubeDataSetToGeoDataFrame(vtk_data)
        gdf.to_file(filename=self._filename)
        #gdf.to_file(filename=self._filename, layer="vtk", driver="GPKG")

        return 1

    def Write(self):
        self.Modified()
        self.Update()


# TODO
#@smproxy.source(name="NCUBEImageOnTopographySource",
#       label="N-Cube Image On Topography Source")

# To add a reader, we can use the following decorators
#   @smproxy.source(name="PythonCSVReader", label="Python-based CSV Reader")
#   @smhint.xml("""<ReaderFactory extensions="csv" file_description="Numpy CSV files" />""")
# or directly use the "@reader" decorator.
@smproxy.reader(name="NCubeLASReader", label="N-Cube LAS Well Log Reader",
                extensions="las", file_description="LAS files")
@smproperty.xml("""<OutputPort name="Header"     index="0" />""")
@smproperty.xml("""<OutputPort name="Curves"     index="1" />""")
class NCubeLASReader(VTKPythonAlgorithmBase):
    """A reader that reads a LAS Well Log file"""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=2)
        self._filename = None
        self._colname = "DEPTH"
        self._x = 0
        self._y = 0
        self._z = 0
        self._az = 0
        self._dip = -90

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="las", file_description="LAS Well Log files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self.Modified()

    @smproperty.doublevector(name="Location", default_values=[0, 0, 0])
    @smdomain.doublerange()
    def SetLocation(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z
        self.Modified()

    @smproperty.doublevector(name="Azimuth", default_values=0)
    @smdomain.doublerange(min=0, max=360)
    def SetAzimuth(self, az):
        self._az = az
        self.Modified()

    @smproperty.doublevector(name="Dip", default_values=-90)
    @smdomain.doublerange(min=-90, max=90)
    def SetDip(self, dip):
        self._dip = dip
        self.Modified()


    def FillOutputPortInformation(self, port, info):
        from vtk import vtkDataObject
        if port == 1:
            info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkPolyData")
        else:
            info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkTable")
        return 1


    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkTable
        #from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtk.util import numpy_support
        from vtk import vtkPolyData, vtkPoints, vtkCellArray, vtkFloatArray, VTK_FLOAT
        import lasio
        import pandas as pd
        import numpy as np
        import math
        import time

        t0 = time.time()
        las = lasio.read(self._filename)
        # DEPTH is index by default
        df_curves = las.df().reset_index()
        headers = []
        for (section, items) in las.sections.items():
            if items is None or items in ('',[]):
                continue
            for item in items:
                headers.append((section,item['mnemonic'],item['unit'],item['value'],item['descr']))
        df_header = pd.DataFrame(headers, columns=('Section','Mnemonic','Unit','Value','Description'))
        vtk_arrays = _NcubeDataFrameToVTKArrays(df_header)
        vtk_table_header = vtkTable()
        for vtk_arr in vtk_arrays:
            vtk_table_header.AddColumn(vtk_arr)

        # https://github.com/mobigroup/gis-snippets/blob/master/ParaView/ProgrammableFilter/vtkMultiblockDataSet.md
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system
        # Spherical coordinates (r, θ, φ) as often used in mathematics:
        # radial distance r, azimuthal angle θ, and polar angle φ.
        theta = 1./2*math.pi - math.pi*self._az/180
        phi = math.pi*(90 - self._dip)/180
        print ("theta",theta,"phi",phi)
        df_curves['dx'] = np.round(df_curves[self._colname]*np.sin(phi)*np.cos(theta),10)
        df_curves['dy'] = np.round(df_curves[self._colname]*np.sin(phi)*np.sin(theta),10)
        df_curves['dz'] = np.round(df_curves[self._colname]*np.cos(phi),10)

        vtk_polyData = vtkPolyData()
        vtk_points = vtkPoints()
        vtk_cells = vtkCellArray()
        vtk_cells.InsertNextCell(len(df_curves))
        for row in df_curves.itertuples(index=False):
            pointId = vtk_points.InsertNextPoint(self._x+row.dx, self._y+row.dy, self._z+row.dz)
            vtk_cells.InsertCellPoint(pointId)
        vtk_polyData.SetPoints(vtk_points)
        vtk_polyData.SetLines(vtk_cells)

        # set of vtk arrays with column names
        vtk_arrays = _NcubeDataFrameToVTKArrays(df_curves)
        for vtk_arr in vtk_arrays:
#            vtk_polyData.GetCellData().AddArray(vtk_arr)
            vtk_polyData.GetPointData().AddArray(vtk_arr)
#        vtk_polyData.GetPointData().SetActiveScalars("DEPTH")

#        print (df_curves['DEPTH'].dtype)
#        vtk_array = numpy_support.numpy_to_vtk(num_array=df_curves['DEPTH'].values, deep=True, array_type=VTK_FLOAT)
#        vtk_array.SetName("DEPTH")
#        vtk_polyData.GetPointData().SetScalars(vtk_array)

        outputHeader = vtkTable.GetData(outInfoVec, 0)
        outputCurves = vtkPolyData.GetData(outInfoVec, 1)

        outputHeader.ShallowCopy(vtk_table_header)
        outputCurves.ShallowCopy(vtk_polyData)

        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1


#@smproxy.filter(name="NCubeTableToLineFilter",
#       label="N-Cube Table To Line Filter")
#@smproperty.input(name="InputTable", port_index=0)
#@smdomain.datatype(dataTypes=["vtkTable"], composite_data_supported=False)
#class NCubeTableToLineFilter(VTKPythonAlgorithmBase):
#    def __init__(self):
#        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, inputType="vtkTable", nOutputPorts=1, outputType="vtkPolyData")
#
#    def RequestData(self, request, inInfoVec, outInfoVec):
#        from vtkmodules.vtkCommonDataModel import vtkTable, vtkPolyData
#        from vtk.util.numpy_support import vtk_to_numpy
#
#        input0 = vtkTable.GetData(inInfoVec[0], 0)
#        output = vtkPolyData.GetData(outInfoVec, 0)
#        # do work
#        print("Pretend work done!")
#        col = input.GetColumnByName('Field 0')
#        print (vtk_to_numpy(col))
#        return 1

