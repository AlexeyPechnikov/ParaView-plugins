# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 

def _NCubeImageOnTopographyToGrid(dem, image):
    from vtk import vtkPoints, vtkStructuredGrid, vtkThreshold, vtkDataObject, VTK_FLOAT, VTK_UNSIGNED_CHAR
    from vtk.util import numpy_support as vn
    import numpy as np

    # mask black and NaN areas
    nanmask = (~np.any(np.isnan(image.values),axis=0)).astype(float)
    # that's correct way
    zeromask = (~np.all(image.values==0,axis=0)).astype(float)
    # that's magic for better borders
    zeromask = (~np.any(image.values==0,axis=0)).astype(float)
    mask = nanmask*zeromask
    mask[mask==0] = np.nan

    xs = dem.x.values
    ys = dem.y.values
    values = mask * dem.values

    colors = np.round(image.values)

    # create raster mask by geometry and for NaNs
    (yy,xx) = np.meshgrid(ys, xs)
    vtk_points = vtkPoints()
    points = np.column_stack((xx.ravel('F'),yy.ravel('F'),values.ravel('C')))
    _points = vn.numpy_to_vtk(points, deep=True)
    vtk_points.SetData(_points)

    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(len(xs), len(ys), 1)
    sgrid.SetPoints(vtk_points)

    array = vn.numpy_to_vtk(values.ravel(), deep=True, array_type=VTK_FLOAT)
    array.SetName("z")
    sgrid.GetPointData().AddArray(array)

    array = vn.numpy_to_vtk(colors.reshape(3,-1).T, deep=True, array_type=VTK_UNSIGNED_CHAR)
    array.SetName("colors")
    sgrid.GetPointData().AddArray(array)

    thresh = vtkThreshold()
    thresh.SetInputData(sgrid)
    thresh.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_POINTS, "z")
    thresh.ThresholdBetween(-1e30, 1e30)
    thresh.Update()

#    return sgrid
    return thresh.GetOutput()

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

        # load the full image raster
        image = xr.open_rasterio(self._imagename)
        image = image.interp_like(dem)
        #dem = dem.interp_like(image)

        vtk_ugrid = _NCubeImageOnTopographyToGrid(dem, image)

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
