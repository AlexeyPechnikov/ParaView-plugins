# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 
from NCube import _NCubeImageOnTopographyToGrid

# load error fix for paraView 5.8.1rc1 Python3
try:
    import xarray
except:
    import sys
    print (sys.exc_info()[0])

#------------------------------------------------------------------------------
# N-Cube Image On Topography Source
#------------------------------------------------------------------------------
@smproxy.source(name="NCubeImageOnTopographySource",
       label="N-Cube Image On Topography Source")
class NCubeImageOnTopographySource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                nInputPorts=0,
                nOutputPorts=1,
                outputType='vtkUnstructuredGrid')
        self._imagename = None
        self._toponame = None
        self._usesealevel = 0
        self._mask_magic = 1


    def RequestData(self, request, inInfo, outInfo):
        from vtk import vtkUnstructuredGrid
        import xarray as xr
        import numpy as np
        import time

        if self._toponame is None or self._imagename is None:
            return 1

        t0 = time.time()

        # load the full topography raster
        dem = xr.open_rasterio(self._toponame).squeeze()
        if dem.values.dtype not in [np.dtype('float16'),np.dtype('float32'),np.dtype('float64'),np.dtype('float128')]:
            dem.values = dem.values.astype("float32")
        dem.values[dem.values == dem.nodatavals[0]] = np.nan
        if self._usesealevel:
            dem.values[dem.values <= 0] = 0

        # load the full image raster
        image = xr.open_rasterio(self._imagename)
        image = image.interp_like(dem)
        #dem = dem.interp_like(image)

        vtk_ugrid = _NCubeImageOnTopographyToGrid(dem, image, self._mask_magic)

        output = vtkUnstructuredGrid.GetData(outInfo, 0)
        output.ShallowCopy(vtk_ugrid)

        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1

    @smproperty.stringvector(name="Image File Name")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["tif", "TIF", "nc"], file_description="GeoTIFF, NetCDF")
    def SetShapeFileName(self, name):
        """Specify filename for the image to read."""
        print ("SetImageFileName", name)
        name = name if name != 'None' else None
        if self._imagename != name:
            self._imagename = name
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

    @smproperty.xml("""
        <IntVectorProperty name="Use Sea Level For Negative Topography"
                       command="SetTopographySeaLevel"
                       number_of_elements="1"
                       default_values="0">
        <BooleanDomain name="bool" />
        <Documentation>
            Use this checkbox to replace negative topography by sea level.
        </Documentation>
        </IntVectorProperty>
    """)
    def SetTopographySeaLevel(self, value):
        print ("TopographySeaLevel", value)
        self._usesealevel = value
        self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Use Magic Image Mask"
                       command="SetUseImageMagicMask"
                       number_of_elements="1"
                       default_values="1">
        <BooleanDomain name="bool" />
        <Documentation>
            Unset this checkbox when you see some missed pixels.
        </Documentation>
        </IntVectorProperty>
    """)
    def SetUseImageMagicMask(self, value):
        print ("SetImageMagicMask", value)
        self._mask_magic = value
        self.Modified()
