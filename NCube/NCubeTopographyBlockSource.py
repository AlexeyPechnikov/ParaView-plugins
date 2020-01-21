# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import * 

from NCube import _NCubeGeoDataFrameLoad, _NCubeTopography, _NCubeTopographyToGrid

# this import is required to fix the issue: https://github.com/Toblerity/Shapely/issues/553
from shapely.geometry import box
import rasterio
import xarray as xr
import numpy as np
import geopandas as gpd
from vtk import vtkCompositeDataSet, vtkMultiBlockDataSet
import time

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
            # TODO: check NODATA in the function
            #dem.values[dem.values == dem.nodatavals[0]] = np.nan
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

