# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

# process [multi]geometry
def _NCubeGeometryToPolyData(geometry, dem=None):
    #from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    import xarray as xr
    import numpy as np

    if geometry is None or geometry.is_empty:
        return

    vtk_points = vtkPoints()
    vtk_cells = vtkCellArray()
    # get part(s) of (multi)geometry
    #if isinstance(geometry, (BaseMultipartGeometry)):
    if geometry.type.startswith('Multi') or geometry.type == 'GeometryCollection':
        geometries = [geom for geom in geometry]
    else:
        geometries = [geometry]
    for geom in geometries:
        # polygon
        #print ("geom.type", geom.type)
        if geom.type == 'Polygon':
            coords = np.asarray(geom.exterior.coords)
        else:
            coords = np.asarray(geom.coords)
        #print ("coords", coords)
        xs = coords[:,0]
        ys = coords[:,1]
        if coords.shape[1] > 2:
            zs = np.array(coords[:,2])
        else:
            zs = np.zeros(len(xs))
        #print (xs)
        # rasterize geometries (lines only, not points)
        # alas, modern scipy or matplotlib don't work in ParaView 5.7 on MacOS
        if dem is not None:
#            print (dem)
            if dem.res and len(xs)>1:
                res = min(dem.res)
                _xs = [xs[:1]]
                _ys = [ys[:1]]
                _zs = [zs[:1]]
                for (x0,y0,z0,x,y,z) in zip(xs[:-1],ys[:-1],zs[:-1],xs[1:],ys[1:],zs[1:]):
                    length = max(abs(x-x0),abs(y-y0))
                    num = round(length/res+0.5)
#                    print ("num",num)
                    if num > 1:
                        _x = np.linspace(x0,x,num)
                        _y = np.linspace(y0,y,num)
                        _z = np.linspace(z0,z,num)
                        _xs.append(_x[1:])
                        _ys.append(_y[1:])
                        _zs.append(_z[1:])
                    else:
                        _xs.append([x])
                        _ys.append([y])
                        _zs.append([z])
                xs = np.concatenate(_xs)
                ys = np.concatenate(_ys)
                zs = np.concatenate(_zs)
            zs += dem.sel(x=xr.DataArray(xs), y=xr.DataArray(ys), method='nearest').values

        #print ("xs", xs)
        mask = np.where(~np.isnan(zs))[0]
        mask2 = np.where(np.diff(mask)!=1)[0]+1
        xs = np.split(xs[mask], mask2)
        ys = np.split(ys[mask], mask2)
        zs = np.split(zs[mask], mask2)
        for (_xs,_ys,_zs) in zip(xs,ys,zs):
            # need to have 2 point or more
            #if len(_xs) <= 1:
            #    continue
            vtk_cells.InsertNextCell(len(_xs))
            for (x,y,z) in zip(_xs,_ys,_zs):
                pointId = vtk_points.InsertNextPoint(x, y, z)
                vtk_cells.InsertCellPoint(pointId)

    # not enougth valid points
    if vtk_points.GetNumberOfPoints() < 1:
        return

    #print ("GetNumberOfPoints", vtk_points.GetNumberOfPoints())
    vtk_polyData = vtkPolyData()
    vtk_polyData.SetPoints(vtk_points)
    #if geometry.type in ['Point','MultiPoint']:
    if geometry.type.endswith('Point'):
        vtk_polyData.SetVerts(vtk_cells)
    else:
        vtk_polyData.SetLines(vtk_cells)
    return vtk_polyData

# process geodataframe and xarray raster
def _NCubeGeometryOnTopography(df, dem):
    from vtk import vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    #import xarray as xr
    import numpy as np

    #print ("_NCUBEGeometryOnTopography start")

    #dem_extent = dem_crs = None
    if dem is not None:
        # TODO: that's better to direct use NODATA values
        if dem.values.dtype not in [np.dtype('float16'),np.dtype('float32'),np.dtype('float64'),np.dtype('float128')]:
            dem.values = dem.values.astype("float32")
        # dask array can't be processed by this way
        dem.values[dem.values == dem.nodatavals[0]] = np.nan
        # NaN border to easy lookup
        dem.values[0,:]  = np.nan
        dem.values[-1,:] = np.nan
        dem.values[:,0]  = np.nan
        dem.values[:,-1] = np.nan

        dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())
        #dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None
        # TODO: FutureWarning: '+init=<authority>:<code>' syntax is deprecated. '<authority>:<code>' is the preferred method.
        dem_crs = dem.crs.replace('+init=','')
        #print (dem.values)
        df = _NCubeGeoDataFrameToTopography(df, dem_extent, dem_crs)

    groups = df.index.unique() ;#[11454:11455]
    #print ("groups",groups)

    # TEST
    #groups = groups[:1]

    # iterate blocks
    vtk_blocks = []
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
#        if hasattr(group, 'encode'):
#            # select only equals
#            _df = df[df.index.str.startswith(group)&df.index.str.endswith(group)&(df.index.str.len()==len(group))].reset_index()
#        else:
        _df = df[df.index == group].reset_index()
        #print (_df.geometry)
        vtk_appendPolyData = vtkAppendPolyData()
        # iterate rows with the same attributes and maybe multiple geometries
        for rowidx,row in _df.iterrows():
            #print ("row", row)
            vtk_polyData = _NCubeGeometryToPolyData(row.geometry, dem)
            if vtk_polyData is None:
                #print ("vtk_polyData is None")
                continue
            vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays(row.to_dict())
            for (vtk_arr, val) in vtk_arrays:
                if val is None:
                    continue
#                for _ in range(vtk_polyData.GetNumberOfCells()):
#                    vtk_arr.InsertNextValue(val)
                if isinstance(val, (tuple)):
#                    if np.any(np.isnan(val)):
#                        continue
                    # add vector
                    for _ in range(vtk_polyData.GetNumberOfCells()):
                        vtk_arr.InsertNextTuple(val)
                    vtk_polyData.GetCellData().AddArray(vtk_arr)
                else:
                    # add scalar
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

        vtk_blocks.append((str(group),vtk_block))

    #print ("_NCUBEGeometryOnTopography end")

    return vtk_blocks

def _NCubeGeoDataFrameToTopography(df, dem_extent, dem_crs=None):
    import geopandas as gpd

    # extract the geometry coordinate system
    df_crs =  df.crs
    print ("df_crs",df_crs,"dem_crs",dem_crs)

    # reproject when the both coordinate systems are defined and these are different
    if df_crs and dem_crs:
        df_extent = gpd.GeoDataFrame([], crs=dem_crs, geometry=[dem_extent])
        print ("df_extent", df_extent.crs, df_extent.geometry)
        extent_reproj = df_extent.to_crs(df_crs)['geometry'][0]
        # if original or reprojected raster extent is valid, use it to crop geometry
        print ("crop geometry", extent_reproj.is_valid,extent_reproj.wkt)
        if extent_reproj.is_valid:
            # geometry intersection to raster extent in geometry coordinate system
            df = df[df.geometry.intersects(extent_reproj)].copy()
            # dangerous operation, see https://github.com/Toblerity/Shapely/issues/553
            df['geometry'] = df.geometry.intersection(extent_reproj)
        return df.to_crs(dem_crs)

    # let's assume the coordinate systems are the same
    if dem_extent is not None:
        df = df[df.geometry.intersects(dem_extent)]
        # wrap issue with 3D geometry intersection by 2D extent
#        if df.geometry[0].has_z:
#            print ("df.geometry[0].has_z")
#        else:
#            df['geometry'] = df.geometry.intersection(dem_extent)

    return df

# Load shapefile or geojson
def _NCubeGeoDataFrameLoad(shapename, shapecol=None, shapeencoding=None):
    import geopandas as gpd

    df = gpd.read_file(shapename, encoding=shapeencoding)
    # very important check
    df = df[df.geometry.notnull()]
    if shapecol is not None:
        df = df.sort_values(shapecol).set_index(shapecol)
    else:
        # to merge all geometries in output
        df.index = len(df)*['None']
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
# we ignore case of scientific notation for numbers
# https://re-thought.com/how-to-suppress-scientific-notation-in-pandas/
def _NCubeGeoDataFrameRowToVTKArrays(items):
    #vtkPolyData, vtkAppendPolyData, vtkPoints, vtkCellArray, 
    from vtk import vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry

    vtk_row = []
    for (key,value) in items.items():
        #print (key,value)
        components = 1
        # define attribute as array
        if isinstance(value, (BaseMultipartGeometry)):
            #print ('BaseMultipartGeometry')
            continue
        elif isinstance(value, (BaseGeometry)):
            #print ('BaseGeometry')
            continue
        elif isinstance(value, (tuple)):
            #print ('vtkFloatArray')
            vtk_arr = vtkFloatArray()
            components = len(value)
#        elif isinstance(value, (int)) or (type(value)==str and value.replace('-','',1).isdigit()):
        elif isinstance(value, (int)) \
                or (type(value)==str and value[0] in ['-','+'] and value[1:].isdigit()) \
                or (type(value)==str and value.isdigit()):
            # ParaView category editor converts strings to numeric when it's possible
            #print('vtkIntArray')
            value = int(value)
            vtk_arr = vtkIntArray()
#        elif isinstance(value, (float)) or (type(value)==str and value.replace('-','',1).replace('.','',1).isdigit()):
        elif isinstance(value, (float)) \
                or (type(value)==str and value[0] in ['-','+'] and value[1:].replace('.','',1).isdigit()) \
                or (type(value)==str and value.replace('.','',1).isdigit()):
            # ParaView category editor converts strings to numeric when it's possible
            #print ('vtkFloatArray')
            value = float(value)
            vtk_arr = vtkFloatArray()
        elif isinstance(value, (bool)):
            #print ('vtkBitArray')
            vtk_arr = vtkBitArray()
        else:
            # some different datatypes could be saved as strings
            value = str(value)
            vtk_arr = vtkStringArray()

        vtk_arr.SetNumberOfComponents(components)
        vtk_arr.SetName(key)
        vtk_row.append((vtk_arr, value))
    return vtk_row

#------------------------------------------------------------------------------
# N-Cube Topography Block Source
#------------------------------------------------------------------------------
def _NCubeTopographyToGrid(dem):
    from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, VTK_FLOAT
    from vtk.util import numpy_support as vn
    import numpy as np
    
    #print (dem)

    xs = dem.x.values
    ys = dem.y.values
    if dem.dims == ('y','x'):
        values = dem.values
    else:
        values = dem.values.T

    # create raster mask by geometry and for NaNs
    #mask = values.ravel('C')==dem.nodatavals[0]
    #print (mask)
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

    # mask NODATA values
    _values = values.ravel('C').astype(float)
    _values[_values == dem.nodatavals[0]] = np.nan
    #array = vn.numpy_to_vtk(values.ravel('C'), deep=True, array_type=VTK_FLOAT)
    array = vn.numpy_to_vtk(_values, deep=True, array_type=VTK_FLOAT)
    array.SetName("z")
    sgrid.GetPointData().AddArray(array)

    thresh = vtkThreshold()
    thresh.SetInputData(sgrid)
    thresh.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, "z")
    thresh.ThresholdBetween(-1e30, 1e30)
    thresh.Update()

#    return sgrid
    return thresh.GetOutput()

# process rasterio dem object plus geodataframe
def _NCubeTopography(dem, df):
    from vtk import vtkAppendFilter
    #vtkPoints, vtkCellArray, vtkStringArray, vtkIntArray, vtkFloatArray, vtkBitArray
    from shapely.geometry.base import BaseGeometry, BaseMultipartGeometry
    from shapely.geometry import box
    import xarray as xr
    import numpy as np
    import rasterio.mask

    dem_crs = dem.crs.to_string() if dem.crs is not None else None
    dem_extent = box(dem.bounds.left,dem.bounds.bottom,dem.bounds.right,dem.bounds.top)
    df = _NCubeGeoDataFrameToTopography(df, dem_extent, dem_crs)

    groups = df.index.unique()
    #print ("groups",groups)

    vtk_blocks = []
    # iterate blocks
    for group in groups:
        #print ("group",group)
        # Python 2 string issue wrapped
#        if hasattr(group, 'encode'):
#            geoms = df[df.index.str.match(group)].geometry
#        else:
        geoms = df[df.index == group].geometry

        vtk_append = vtkAppendFilter()
        # rasterize geomemetry
        try:
            # use only 1st band from the raster
            out_image, out_transform = rasterio.mask.mask(dem, geoms, indexes=1, crop=True, filled=True)
        except:
            # geometry outside of the raster, etc.
            print (group, "exception")
            continue
        (xs, _) = rasterio.transform.xy(out_transform, out_image.shape[1]*[0], range(out_image.shape[1]), offset='center')
        (_, ys) = rasterio.transform.xy(out_transform, range(out_image.shape[0]), out_image.shape[0]*[0], offset='center')
        da = xr.DataArray(out_image, coords={'x': xs, 'y': ys}, dims=('y', 'x'))
        da.attrs['nodatavals'] = (dem.nodata,)
        # process rasterized geometry
        vtk_ugrid = _NCubeTopographyToGrid(da)
        if vtk_ugrid is None:
            continue
        if df.index is not None:
            vtk_arrays = _NCubeGeoDataFrameRowToVTKArrays({df.index.name:group})
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
        vtk_blocks.append((str(group),vtk_block))

    print ("_NCubeTopography end")

    return vtk_blocks

#------------------------------------------------------------------------------
# N-Cube Image On Topography Source
#------------------------------------------------------------------------------
def _NCubeImageOnTopographyToGrid(dem, image, mask_magic=False):
    from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, VTK_FLOAT, VTK_UNSIGNED_CHAR
    from vtk.util import numpy_support as vn
    import numpy as np

    # mask NaN areas
    nanmask = (~np.any(np.isnan(image.values),axis=0)).astype(float)
    # mask single channel zeroes if needed
    if mask_magic:
        # that's more correct way, only black pixels ignored
        #zeromask = (~np.all(image.values==0,axis=0)).astype(float)
        # that's magic for better borders
        zeromask = (~np.any(image.values==0,axis=0)).astype(float)
        mask = nanmask*zeromask
    else:
        mask = nanmask
    mask[mask==0] = np.nan

    xs = dem.x.values
    ys = dem.y.values
    values = mask * dem.values

    # create raster mask by geometry and for NaNs
    (yy,xx) = np.meshgrid(ys, xs)
    vtk_points = vtkPoints()
    points = np.column_stack((xx.ravel('F'),yy.ravel('F'),values.ravel('C')))
    _points = vn.numpy_to_vtk(points, deep=True)
    vtk_points.SetData(_points)

    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(len(xs), len(ys), 1)
    sgrid.SetPoints(vtk_points)

    # convert datatype to Float32
    array = vn.numpy_to_vtk(values.astype(np.float32).ravel(), deep=True, array_type=VTK_FLOAT)
    array.SetName("z")
    sgrid.GetPointData().AddArray(array)

    # iterate 2D only coordinates
    for coord in dem.coords:
        if len(dem.coords[coord].shape) != 2:
            continue
        print (len(dem.coords[coord].shape))
        # convert datatype to Float32
        array = vn.numpy_to_vtk(dem.coords[coord].values.astype(np.float32).ravel(), deep=True, array_type=VTK_FLOAT)
        array.SetName(coord)
        sgrid.GetPointData().AddArray(array)

    bands = image.band.shape[0]
    print ("bands", bands)
    if bands == 3:
        # RGB
        colors = np.round(image.values)
        array = vn.numpy_to_vtk(colors.reshape(3,-1).T, deep=True, array_type=VTK_UNSIGNED_CHAR)
        array.SetName("colors")
        sgrid.GetPointData().AddArray(array)
    elif bands == 1:
        arr = image.values
        array = vn.numpy_to_vtk(arr.reshape(1,-1).T, deep=True, array_type=VTK_FLOAT)
        array.SetName("band")
        sgrid.GetPointData().AddArray(array)
    else:
        print ("Unsupported bands count (should be 1 or 3)", bands)

    thresh = vtkThreshold()
    thresh.SetInputData(sgrid)
    thresh.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, "z")
    thresh.ThresholdBetween(-1e30, 1e30)
    thresh.Update()

#    return sgrid
    return thresh.GetOutput()

