# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import * 

from NCube import _NCubeGeoDataFrameLoad, _NCubeGeoDataFrameToTopography, _NCubeGeoDataFrameRowToVTKArrays

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
        if hasattr(group, 'encode'):
            geoms = df[df.index.str.match(group)].geometry
        else:
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
        # this import is required to fix the issue: https://github.com/Toblerity/Shapely/issues/553
        from shapely.geometry import box
        import rasterio
        import xarray as xr
        import numpy as np
        import geopandas as gpd
        from vtk import vtkCompositeDataSet, vtkMultiBlockDataSet
        import time

        if self._toponame is None:
            return 1

        t0 = time.time()

        # load geometries
        if self._shapename is not None:
            df = _NCubeGeoDataFrameLoad(self._shapename, self._shapecol, self._shapeencoding)
        else:
            df = None

        if df is None:
            # process the full topography raster
            dem = xr.open_rasterio(self._toponame).squeeze()
            vtk_ugrid = _NCubeTopographyToGrid(dem)
            #print ("vtk_ugrid",vtk_ugrid)
            vtk_blocks = [(str('None'),vtk_ugrid)]
        else:
            # open raster
            dem = rasterio.open(self._toponame)
            # process shapefile
            vtk_blocks = _NCubeTopography(dem, df)
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
        import geopandas as gpd
        if self._shapename is None:
            return []
        # Load shapefile
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

