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
def _NCubeGeoDataFrameLoad(shapename, shapecol=None, shapeencoding=None):
    import geopandas as gpd
    if shapename is None:
        return
    df = gpd.read_file(shapename, encoding=shapeencoding)
    # clean undefined geometries
    if (len(df)) == 0:
        return []
    # pre-generate sindex on df1 if it doesn't already exist
    df.sindex
    if shapecol is not None:
        return df.sort_values(shapecol).set_index(shapecol)
    #print ("shapecol",shapecol)
    return df

# load DEM
def _NCubeRasterLoad(rastername):
    import xarray as xr
    import numpy as np
    if rastername is None:
        return
    raster = xr.open_rasterio(rastername).squeeze()
    raster.values[raster.values == raster.nodatavals[0]] = np.nan
    #raster.values[raster.values == raster.nodatavals[0]] = 0

    # NaN border to easy lookup
    raster.values[0,:]  = np.nan
    raster.values[-1,:] = np.nan
    raster.values[:,0]  = np.nan
    raster.values[:,-1] = np.nan

    return raster

def _NCubeGeoDataFrameToRaster(df, extent, dem_crs):
    #from shapely.geometry import Polygon
    from functools import partial
    import pyproj
    from shapely.ops import transform

    if df is None:
        return

    # reproject when the both coordinate systems are defined and these are different
    df_crs = str(df.crs['init']) if df.crs != {} else None
    if df_crs and dem_crs:
        print ("df_crs",df_crs,"dem_crs",dem_crs)
        # for Python2
        project = partial(
            pyproj.transform,
            pyproj.Proj(init=df_crs), # source coordinate system
            pyproj.Proj(init=dem_crs) # destination coordinate system
        )
        extent_reproj = transform(project, extent)
        # if original or reprojected raster extent is valid, use it to crop geometry
        if extent_reproj.is_valid:
            # geometry intersection to raster extent in geometry coordinate system
            df = df[df.geometry.intersects(extent_reproj)]
            df['geometry'] = df.geometry.intersection(extent_reproj)

        # reproject [cropped] geometry to original raster coordinates if needed
        df['geometry'] = [transform(project, geom) if geom.is_valid else None for geom in df['geometry']]
        # fix for broken [reprojected] geometries
        df = df[df.geometry.is_valid==True].copy()

    return df


# list of list of VtkArray's
def _NCubeGeoDataFrameRowToVTKArrays(row):
    #vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, 
    from vtk import vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry

    vtk_row = []
    for (key,value) in row.to_dict().items():
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


def _NCubeGeometryToPolyData(geometry, dem=None):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    import xarray as xr
    import numpy as np

    if isinstance(geometry, (BaseMultipartGeometry)):
    #if row.geometry.geometryType()[:5] == 'Multi':
        geoms = geometry.geoms
    else:
        geoms = [geometry]

    vtk_points = vtkPoints()
    vtk_cells = vtkCellArray()
    # iterate parts of (multi)geometry
    for geom in geoms:
        if geom.type == 'Polygon':
            # use exterior coordinates only
            coords = geom.exterior.coords
        else:
            coords = geom.coords

        xs = np.array(coords.xy[0])
        ys = np.array(coords.xy[1])
        if len(coords.xy) > 2:
            zs = np.array(coords.xy[2])
        else:
            zs = np.array([0]*len(xs))
        # find nearest raster values for every geometry point
        if dem is not None:
            zs = dem.sel(x=xr.DataArray(xs), y=xr.DataArray(ys), method='nearest').values
            # code below should be faster for large geometries
#            yy = dem.y.sel(y=xr.DataArray(ys), method='nearest').values
#            xx = dem.x.sel(x=xr.DataArray(xs), method='nearest').values
#            zs = dem.sel(x=xr.DataArray(xx), y=xr.DataArray(yy)).values
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
#    if geom.type == 'Point':
#        vtk_polyData.SetVerts(vtk_cells)
#    else:
    vtk_polyData.SetLines(vtk_cells)

    return vtk_polyData


def _NCubeGeometryOnTopography(shapename, toponame, shapecol, shapeencoding):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    import numpy as np

    #print ("_NCUBEGeometryOnTopography start")

    # Load shapefile
    if shapename is None:
        return

    df = _NCubeGeoDataFrameLoad(shapename, shapecol, shapeencoding)
    if df is None or len(df) == 0:
        return

    # load DEM
    dem = None
    if toponame is not None:
        dem = _NCubeRasterLoad(toponame)
        if dem is None:
            return
    dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())
    dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None

    # crop geometry and topography together
    df = _NCubeGeoDataFrameToRaster(df, dem_extent, dem_crs)

    groups = df.index.unique()
    #print ("groups",groups)

    # iterate blocks
    vtk_blocks = []
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
        if hasattr(group, 'encode'):
            _df = df[df.index.str.match(group)].reset_index()
        else:
            _df = df[df.index == group].reset_index()
        vtk_appendPolyData = vtkAppendPolyData()
        # iterate rows
        for rowidx,row in _df.iterrows():
            vtk_polyData = _NCubeGeometryToPolyData(row.geometry, dem)
            if vtk_polyData is None:
                continue
            vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays(row)
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
def _NCubeTopographyToPolyData(dem, geometry=None):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkCleanPolyData, vtkPoints, vtkCellArray, vtkTriangle, vtkQuad, vtkFloatArray
    from shapely.geometry import Point
    import numpy as np

    xs = dem.x.values
    ys = dem.y.values
    if dem.dims == ('y','x'):
        values = dem.values
    else:
        values = dem.values.T

    # create raster mask by geometry and for NaNs
    (yy,xx) = np.meshgrid(ys, xs)
    mask = ~np.isnan(values)
    if geometry is not None:
        (xmin,ymin,xmax,ymax) = geometry.bounds
        mask = [Point(x,y).intersects(geometry) if m and (x>=xmin and x<=xmax and y>=ymin and y<=ymax) else False for (x,y,m) in zip(xx.ravel('F'), yy.ravel('F'), mask.flatten())]
        mask = np.array(mask).reshape(len(ys),len(xs))
    # nothing to do: actually we need 4+ points to build 1 cell
    if mask.sum() == 0:
        return

    # Define points and triangles for mesh
    vtk_points = vtkPoints()
    vtk_cells = vtkCellArray()

    # Build the meshgrid manually
    count = 0
    # iterate array
    for j in range(0,len(ys)-2):
        for i in range(0,len(xs)-2):
            # check area
            if not mask[j,i] or not mask[j+1,i] or not mask[j,i+1] or not mask[j+1,i+1]:
                continue

            vtk_points.InsertNextPoint(xs[i],   ys[j],   values[j,i])
            vtk_points.InsertNextPoint(xs[i],   ys[j+1], values[j+1,i])
            vtk_points.InsertNextPoint(xs[i+1], ys[j+1], values[j+1,i+1])
            vtk_points.InsertNextPoint(xs[i+1], ys[j],   values[j,i+1])

            quad = vtkQuad()
            quad.GetPointIds().SetId(0, count)
            quad.GetPointIds().SetId(1, count + 1)
            quad.GetPointIds().SetId(2, count + 2)
            quad.GetPointIds().SetId(3, count + 3)

            vtk_cells.InsertNextCell(quad)

            count += 4

    # nothing to do
    if count == 0:
        return

    # Create a polydata object
    trianglePolyData = vtkPolyData()

    # Add the geometry and topology to the polydata
    trianglePolyData.SetPoints(vtk_points)
    #trianglePolyData.GetPointData().SetScalars(colors)
    trianglePolyData.SetPolys(vtk_cells)

    # Clean the polydata so that the edges are shared
    cleanPolyData = vtkCleanPolyData()
    cleanPolyData.SetInputData(trianglePolyData)
    cleanPolyData.Update()

    output = cleanPolyData.GetOutput()
    array = vtkFloatArray()
    array.SetName("z")
    for i in range(0, output.GetNumberOfPoints()):
        array.InsertNextValue(output.GetPoint(i)[2])
    output.GetPointData().SetScalars(array)

    return output


def _NCubeTopography(shapename, toponame, shapecol, shapeencoding):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    import numpy as np

    print ("_NCubeTopography start")

    # Load shapefile
    if toponame is None:
        return

    df = None
    if shapename is not None:
        df = _NCubeGeoDataFrameLoad(shapename, shapecol, shapeencoding)
        if df is None or len(df) == 0:
            return

    # load DEM
    dem = _NCubeRasterLoad(toponame)
    if dem is None:
        return
    dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None
    dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())

    # process the full topography raster
    if df is None:
        vtk_polyData = _NCubeTopographyToPolyData(dem, geometry=None)
        return [(_str('None'),vtk_polyData)]

    # crop geometry and topography together
    df = _NCubeGeoDataFrameToRaster(df, dem_extent, dem_crs)
    #print ("df", list(df.index))

    groups = df.index.unique()
    #print ("groups",groups)

    # iterate blocks
    vtk_blocks = []
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
        if hasattr(group, 'encode'):
            _df = df[df.index.str.match(group)].reset_index()
        else:
            _df = df[df.index == group].reset_index()
        if shapecol is not None:
            _df = _df[[shapecol,'geometry']]
        vtk_appendPolyData = vtkAppendPolyData()
        # iterate rows
        for rowidx,row in _df.iterrows():
            # processing closed areas only
            if not row.geometry.type in ['Polygon', 'MultiPolygon']:
                continue
            vtk_polyData = _NCubeTopographyToPolyData(dem, row.geometry)
            if vtk_polyData is None:
                continue
            if shapecol is not None:
                vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays(row)
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

    print ("_NCubeTopography end")

    return vtk_blocks


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
        if vtk_blocks == []:
            return
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
        from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
        import time

        if self._toponame is None:
            return 1

        t0 = time.time()
        vtk_blocks = _NCubeTopography(self._shapename, self._toponame, self._shapecol, self._shapeencoding)
        if vtk_blocks == []:
            return
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

# TODO
#@smproxy.source(name="NCUBEImageOnTopographySource",
#       label="N-Cube Image On Topography Source")

