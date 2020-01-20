# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 

from NCube import _NCubeDataSetToGeoDataFrame

from vtkmodules.vtkCommonDataModel import vtkDataSet

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
