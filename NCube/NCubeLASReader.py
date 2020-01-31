# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import * 

from NCube import _NcubeDataFrameToVTKArrays

@smproxy.reader(name="NCubeLASReader", label="N-Cube LAS Well Log Reader",
                extensions="las", file_description="LAS files")
@smproperty.xml("""<OutputPort name="Header"     index="0" />""")
@smproperty.xml("""<OutputPort name="Curves"     index="1" />""")
class NCubeLASReader(VTKPythonAlgorithmBase):
    """A reader that reads a LAS Well Log file"""
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=2)
        self._filename = None
        self._x = 0
        self._y = 0
        self._z = 0
        self._az = 0
        self._dip = -90

    def FillOutputPortInformation(self, port, info):
        from vtk import vtkDataObject
        if port == 1:
            info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkPolyData")
        else:
            info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkTable")
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkTable
        from vtk.util import numpy_support
        from vtk import vtkPolyData, vtkPoints, vtkCellArray, vtkFloatArray, VTK_FLOAT
        import pandas as pd
        import numpy as np
        import lasio
        import math
        import time

        t0 = time.time()
        las = lasio.read(self._filename)

        # DEPTH is index
        df_curves = las.df()
        headers = []
        for (section, items) in las.sections.items():
            if items is None or items in ('',[]):
                continue
            if isinstance(items, (str,unicode)):
                headers.append((section,'','',items,''))
            elif isinstance(items, (list)):
                for item in items:
                    headers.append((section,item['mnemonic'],item['unit'],item['value'],item['descr']))
            else:
                print ("Unknown LAS header section type", type(items), iyems)
        df_header = pd.DataFrame(headers, columns=('Section','Mnemonic','Unit','Value','Description'))
        vtk_arrays = _NcubeDataFrameToVTKArrays(df_header)
        vtk_table_header = vtkTable()
        for vtk_arr in vtk_arrays:
            vtk_table_header.AddColumn(vtk_arr)

        outputHeader = vtkTable.GetData(outInfoVec, 0)
        outputHeader.ShallowCopy(vtk_table_header)

        # define scale factor by depth units
        unit = las.curves[0]['unit']
        if unit in ['FT','ft']:
            scale = 0.3048
        elif unit in ['M','m']:
            scale = 1.0
        else:
            # we can use additional ParaView filter to fix it later
            scale = 1.0
            print ("Unknown LAS header unit", unit)
        # set of vtk arrays with column names
        vtk_arrays = _NcubeDataFrameToVTKArrays(df_curves)
        # https://github.com/mobigroup/gis-snippets/blob/master/ParaView/ProgrammableFilter/vtkMultiblockDataSet.md
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system
        # Spherical coordinates (r, θ, φ) as often used in mathematics:
        # radial distance r, azimuthal angle θ, and polar angle φ.
        theta = 1./2*math.pi - math.pi*self._az/180
        phi = math.pi*(90 - self._dip)/180
        #print ("theta",theta,"phi",phi)
        df_curves['dx'] = np.round(scale*df_curves.index*np.sin(phi)*np.cos(theta),10)
        df_curves['dy'] = np.round(scale*df_curves.index*np.sin(phi)*np.sin(theta),10)
        df_curves['dz'] = np.round(scale*df_curves.index*np.cos(phi),10)

        vtk_polyData = vtkPolyData()
        vtk_points = vtkPoints()
        vtk_cells = vtkCellArray()
        vtk_cells.InsertNextCell(len(df_curves))
        for row in df_curves.itertuples(index=False):
            pointId = vtk_points.InsertNextPoint(self._x+row.dx, self._y+row.dy, self._z+row.dz)
            vtk_cells.InsertCellPoint(pointId)
        vtk_polyData.SetPoints(vtk_points)
        vtk_polyData.SetLines(vtk_cells)

        for vtk_arr in vtk_arrays:
#            vtk_polyData.GetCellData().AddArray(vtk_arr)
            vtk_polyData.GetPointData().AddArray(vtk_arr)
#        vtk_polyData.GetPointData().SetActiveScalars("DEPTH")

#        print (df_curves['DEPTH'].dtype)
#        vtk_array = numpy_support.numpy_to_vtk(num_array=df_curves['DEPTH'].values, deep=True, array_type=VTK_FLOAT)
#        vtk_array.SetName("DEPTH")
#        vtk_polyData.GetPointData().SetScalars(vtk_array)

        outputCurves = vtkPolyData.GetData(outInfoVec, 1)
        outputCurves.ShallowCopy(vtk_polyData)

        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1

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

