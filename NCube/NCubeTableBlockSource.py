# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import * 

from NCube import _NCubeGeometryOnTopography

#------------------------------------------------------------------------------
# N-Cube Topography Block Source
#
# Create 3D points when TOP X, Y, Z Fields
#
# OR
#
# Create 3D lines when TOP X,Y,Z Fields defined and some of Bottom X,Y,Z Fields defined
#
#------------------------------------------------------------------------------

@smproxy.source(name="NCubeTableBlockSource",
       label="N-Cube Table Block Source")
class NCubeTableBlockSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                nInputPorts=0,
                nOutputPorts=1,
                outputType='vtkMultiBlockDataSet')
        self._scaleunits = {'m': 1, 'km': 1000, 'ft': 0.3048, 'kft': 304.8}

        self._tablename = None ;# required
        self._tablecol = None
        self._tableencoding = None
        self._epsg = None

        self._xcol = None ;# required
        self._ycol = None ;# required
        self._zcol = None ;# optional - topography could be used instead

        # if it's defined, create segment, overwise create point
        self._xcol2 = None
        self._ycol2 = None
        self._zcol2 = None


    def RequestData(self, request, inInfo, outInfo):
        import xarray as xr
        import numpy as np
        import geopandas as gpd
        import pandas as pd
        from shapely.geometry import box, Point, LineString
        from shapely.affinity import translate
        from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
        import math
        import time

        # check mandatory fields
        if self._tablename is None or self._xcol is None or self._ycol is None or self._zcol is None:
            return 1

        t0 = time.time()
        #df = pd.read_csv(self._tablename, low_memory=False)

        df = pd.read_csv(self._tablename, encoding=self._tableencoding, low_memory=False)
        print ("len(df)",len(df))
        if (len(df)) == 0:
            return

        # calculate coordinates
        xs = df[self._xcol]
        ys = df[self._ycol]
        zs = df[self._zcol]
        # create point
        geom = gpd.GeoSeries(map(Point, zip(xs, ys, zs)))

        # create line segment when possible
        if self._xcol2 is not None or self._ycol2 is not None or self._zcol2 is not None:
            # we can define only one (changed) column for 2nd point
            self._xcol2 = self._xcol2 if self._xcol2 is not None else self._xcol
            self._ycol2 = self._ycol2 if self._ycol2 is not None else self._ycol
            self._zcol2 = self._zcol2 if self._zcol2 is not None else self._zcol
            xs2 = df[self._xcol2]
            ys2 = df[self._ycol2]
            zs2 = df[self._zcol2]
            # create line segment
            geom2= gpd.GeoSeries(map(Point, zip(xs2, ys2, zs2)))
            geom = gpd.GeoSeries(map(LineString, zip(geom,geom2)))

        # create geodataframe to build output
        df = gpd.GeoDataFrame(df, geometry=geom)

        # group to multiblock
        if self._tablecol is not None:
            df = df.sort_values(self._tablecol).set_index(self._tablecol)

        vtk_blocks = _NCubeGeometryOnTopography(df, None)
        if len(vtk_blocks)>0:
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


    # possible coordinate columns
    @smproperty.stringvector(name="NumericColumns", information_only="1")
    def GetNumericColumns(self):
        import pandas as pd
        print ("GetNumericColumns")
        if self._tablename is None:
            return []
        # Load file
        df = pd.read_csv(self._tablename, nrows=1)
        cols = sorted(df._get_numeric_data().columns)
        #print ("GetColumns", cols)
        return ['None'] + list(map(str,cols))

    # Multiblock labels
    @smproperty.stringvector(name="Columns", information_only="1")
    def GetColumns(self):
        import pandas as pd
        if self._tablename is None:
            return []
        # Load file
        df = pd.read_csv(self._tablename, nrows=0)
        cols = sorted(df.columns.values)
        return ['None'] + list(map(str,cols))


    @smproperty.stringvector(name="TableFile")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["csv"], file_description="CSV")
    def SetShapeFileName(self, name):
        """Specify filename for the table to read."""
        print ("SetCSVFileName", name)
        name = name if name != 'None' else None
        if self._tablename != name:
            self._tablename = name
            self.Modified()

    @smproperty.stringvector(name="Table Group By", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Columns" function="GetColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetTableGroupBy(self, col):
        col = col if col != 'None' else None
        self._tablecol = col
        print("SetTableGroupBy", col)
        self.Modified()

    @smproperty.stringvector(name="Top Easting", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetTopXColumn(self, col):
        col = col if col != 'None' else None
        self._xcol = col
        print("SetXColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Top Northing", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetTopYColumn(self, col):
        col = col if col != 'None' else None
        self._ycol = col
        print("SetYColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Top Elevation", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetTopZColumn(self, col):
        col = col if col != 'None' else None
        self._zcol = col
        print("SetZColumn", col)
        self.Modified()


    @smproperty.stringvector(name="Bottom Easting", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetBottomXColumn(self, col):
        col = col if col != 'None' else None
        self._xcol2 = col
        print("SetXColumn2", col)
        self.Modified()

    @smproperty.stringvector(name="Bottom Northing", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetBottomYColumn(self, col):
        col = col if col != 'None' else None
        self._ycol2 = col
        print("SetYColumn2", col)
        self.Modified()

    @smproperty.stringvector(name="Bottom Elevation", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetBottomZColumn(self, col):
        col = col if col != 'None' else None
        self._zcol2 = col
        print("SetZColumn2", col)
        self.Modified()

