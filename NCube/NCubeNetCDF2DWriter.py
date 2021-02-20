# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 
from vtkmodules.vtkCommonDataModel import vtkDataSet

def _NCubeDataSetToNetCDF2D(vtk_data):
    import pandas as pd
    import numpy as np
    from vtk.util import numpy_support

    # save all points with attributes
    coords = vtk_data.GetPoints().GetData()
    ncoords = numpy_support.vtk_to_numpy(coords)
    df = pd.DataFrame(data=ncoords,columns=('x','y','z'))

    for idx in range(vtk_data.GetPointData().GetNumberOfArrays()):
        col = vtk_data.GetPointData().GetArrayName(idx)
        values = vtk_data.GetPointData().GetArray(idx)
        nvalues = numpy_support.vtk_to_numpy(values)
        #nvalues[np.abs(nvalues)>1e37] = np.nan
        df[col] = nvalues

    ds = df.groupby(['z','y','x']).mean().to_xarray().squeeze()
    #print (ds)

    return ds

#------------------------------------------------------------------------------
# N-Cube DataSet to Shapefile Writer
#------------------------------------------------------------------------------
@smproxy.writer(extensions="nc", file_description="NetCDF 2D", support_reload=False)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class NCubeNetCDF2DWriter(VTKPythonAlgorithmBase):
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
        print ("NCubeNetCDF2DWriter")
        vtk_data = vtkDataSet.GetData(inInfoVec[0], 0)
        ds = _NCubeDataSetToNetCDF2D(vtk_data)
        ds.to_netcdf(self._filename)
        return 1

    def Write(self):
        self.Modified()
        self.Update()
