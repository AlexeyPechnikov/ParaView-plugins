# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 
from vtkmodules.vtkCommonDataModel import vtkDataSet

#from NCube import _NCubeDataSetToGeoDataFrame
def _NCubeDataSetToGeoDataFrame(vtk_data, extent_only=False):
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point, box
    from vtk.util import numpy_support

    # save all the points or only the boundaries
    if extent_only or vtk_data.GetNumberOfPoints() == 0:
        # save only boundaries
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
        vtk_data = vtkDataSet.GetData(inInfoVec[0], 0)
        gdf = _NCubeDataSetToGeoDataFrame(vtk_data)
        gdf.to_file(filename=self._filename)
        #gdf.to_file(filename=self._filename, layer="vtk", driver="GPKG")

        return 1

    def Write(self):
        self.Modified()
        self.Update()

@smproxy.writer(extensions="shp", file_description="ESRI Shapefile - Extent only", support_reload=False)
#@smproxy.writer(extensions="gpkg", file_description="GeoPackage", support_reload=False)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class NCubeShapefileWriter2(VTKPythonAlgorithmBase):
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
        vtk_data = vtkDataSet.GetData(inInfoVec[0], 0)
        gdf = _NCubeDataSetToGeoDataFrame(vtk_data, extent_only=True)
        gdf.to_file(filename=self._filename)
        #gdf.to_file(filename=self._filename, layer="vtk", driver="GPKG")

        return 1

    def Write(self):
        self.Modified()
        self.Update()

