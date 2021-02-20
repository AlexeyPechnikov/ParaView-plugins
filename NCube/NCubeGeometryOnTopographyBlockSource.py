# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import *
from NCube import _NCubeGeoDataFrameLoad, _NCubeGeometryOnTopography

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
        self._ignorez = False


    def RequestData(self, request, inInfo, outInfo):
        import xarray as xr
        import numpy as np
        import geopandas as gpd
        from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
        from shapely.ops import transform
        import time

        if self._shapename is None:
            return 1

        df = _NCubeGeoDataFrameLoad(self._shapename, self._shapecol, self._shapeencoding)
        if df is None:
            return
        if self._ignorez:
            df['geometry'] = df.geometry.apply(lambda geom: transform(lambda x, y, z=None: (x, y), geom))

        # load DEM
        dem = None
        if self._toponame is not None:
            #dem = xr.open_rasterio(toponame, chunks=10000000).squeeze()
            dem = xr.open_rasterio(self._toponame).squeeze()

        t0 = time.time()
        vtk_blocks = _NCubeGeometryOnTopography(df, dem)
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

    @smproperty.stringvector(name="ShapefileName")
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

    @smproperty.xml("""
        <IntVectorProperty name="ZIgnore"
                       command="SetShapeFileZIgnore"
                       number_of_elements="1"
                       default_values="0">
        <BooleanDomain name="bool" />
        <Documentation>
            Ignore Z coordinate in 3D shapefiles.
        </Documentation>
        </IntVectorProperty>
    """)
    def SetShapeFileZIgnore(self, value):
        print ("SetShapeFileZIgnore", value)
        self._ignorez = value
        self.Modified()
