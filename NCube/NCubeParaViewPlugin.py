# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

#from NCube import *
from NCubeGeometryOnTopographyBlockSource import *
from NCubeLASReader import *
from NCubeShapefileWriter import *
from NCubeTableOnTopographyBlockSource import *
from NCubeTopographyBlockSource import *

#from paraview.util.vtkAlgorithm import *
#from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase


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
