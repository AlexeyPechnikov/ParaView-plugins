# N-Cube ParaView plugin for 3D/4D GIS Data Visualization

<img src="https://github.com/mobigroup/gis-snippets/blob/master/ParaView/ProgrammableFilter/ParaView_ProgrammableFilter_reproject.jpg" width="400">

## About the project

[N-Cube ParaView plugin](NCube/NCubeParaViewPlugin.py) is [MIT-licensed](LICENSE) set of ParaView data readers/sources/filters/writers to load, process, and write well-known GIS datasets like to shapefiles, GeoTIFF rasters, etc.

We use this plugin on MacOS only although it's cross platform and can be used with Linux and Windows ParaView builds too. 

## How to Use ParaView 5.9.0 With Python 3.8 on MacOS

Fortunately, ParaView 5.9 on MacOS supports Python 3.8 and we can drop Python 2.7 support now! Today protected MacOS applications (signed and hardened) are not able to load 3rd party binary libraries (required for the plugin) and we need to use special ParaView build to load them: [ParaView-5.9.0-MPI-OSX10.13-Python3.8-64bit.unsigned.dmg](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.9&type=binary&os=macOS&downloadFile=ParaView-5.9.0-MPI-OSX10.13-Python3.8-64bit.unsigned.dmg)

Note: while I tried to sign 3rd party binaries compiled by GCC compiler by self-issued certificate it doesn't work (missed Apple compiler specific LC_VERSION_MIN) and we need to use the special ParaView build. Let me know if you know how to resolve the issue!

To extend your ParaView by 3rd party Python libraries and ParaView Python plugins follow instructions below:

 * [One-time step] Install Python 3.8 using [Homebrew](https://brew.sh):
 ```
 brew install python@3.8
 ```

 * [One-time step] Install required Python modules:
```
/usr/local/opt/python@3.8/bin/python3 -m pip install --upgrade \
  numpy xarray pandas geopandas shapely vtk rasterio lasio
```

* [One-time step] Create ParaView wrapper script pv5.9 with this content:
```
#!/bin/sh
PYTHONPATH=/usr/local/lib/python3.8/site-packages
GDAL_DATA=/usr/local/lib/python3.8/site-packages/rasterio/gdal_data
APPPATH=/Applications/ParaView-5.9.0.app
PYTHONPATH="$PYTHONPATH" GDAL_DATA="$GDAL_DATA" "$APPPATH/Contents/MacOS/paraview"
```

 * Launch ParaView from MacOS Terminal by command below and save the Terminal window open while you use ParaView:
 ```
 pv5.9
 ```
Note: you are able to see debug messages for developers in the Terminal Window.

 * [One-time step] For the first launch download the plugin code or clone the repository and load from ParaView("Tools" -> "Manage Plugins" -> "Load New ...") this file [N-Cube ParaView plugin Python source file](NCube/NCubeParaViewPlugin.py) placed anywhere on your computer. Setup 'Auto Load' checkbox if you need to use the plugin later otherwise it should be loaded at next ParaView launches again.

<img src="screenshots/NcubePlugin.jpg" width="400">

## How to Use old ParaView versions With Python 2.7

If you still need to have Python 2.7 compatible build see this link: [Version compatiable with Python 2.7 and Python 3.7+](https://github.com/mobigroup/ParaView-plugins/releases/tag/python2.7) Use the same instructions as above for Python 3.8.

Note: there are some non-latin symbols issues with Python 2.7. Maybe you need to remove non-latin fields in your data files and rename them to use ASCII symbols only.

### N-Cube ParaView Plugin Readers:

**N-Cube LAS Well Log Reader** - Read Well Log versions 1.2 and 2.0 of the LAS file specification.

### N-Cube ParaView Plugin Writers:

**ESRI Shapefile** - Use ParaView menu Save -> Save Data -> ESRI Shapefile(\*.shp) to save geometry as ESRI Shapefile (Point). By performance reasons in case when the geometry includes more than 1M points only bounding box will be saved.

### N-Cube ParaView Plugin Sources:

**N-Cube Image On Topography Source** - data source for georeferenced image. Optionally GeoTIFF or NetCDF topography (DEM) raster using to define Z coordinates. Optional "Use Sea Level for Negative Topography" parameter allows to replace negative topography values by zeroes.

**N-Cube Reproject Filter** - TODO.

**N-Cube Shapefile On Topography Block Source** - data source for 2D/3D Shapefile or GeoJSON. Optionally GeoTIFF or NetCDF topography (DEM) raster using to define Z coordinates. Optional "Group by Field" parameter produces a set of layers from geometries grouped by the field. When the "Group by Field" parameter is not defined ("None") all the geometries separated as the layers. All shapefile fields presented in the ouput.

**N-Cube Table Block Source** - data source for CSV table data to produce geometries like to wells. 

**N-Cube Table on Topography Block Source** - data source for CSV table data to produce geometries on optional topography. See for details [NCubeTableOnTopographyBlockSource](https://github.com/mobigroup/ParaView-plugins/tree/master/NCube/NCubeTableOnTopographyBlockSource)

**N-Cube Topography Block Source** - data source for GeoTIFF or NetCDF topography (DEM) raster visualization as 3D surface. Only "Group by Field" parameter presented in the ouput.

Geometry coordinate system re-projecting to raster coordinate system and the geometry cropping to the raster extent when the topography file using and the both coordinate systems are defined.

## Project Goals

For our 3D geological modeling we need a good 3D visualization and data processing tools. We tested many commercial and Open Source 3D visualization packages and programming libraries and for our needs Open Source ParaView software and it's core VTK library are the best one. [ParaView](https://www.paraview.org/) is the great 3D visualization and processing tool but without GIS data support. To fix the lack, we created ParaView geospatial plugins for our internal usage.

See also our code snippets repository [ParaView Programmable Source and Programmable Filter examples](https://github.com/mobigroup/gis-snippets/tree/master/ParaView) and [A brief explanation of the 3D Density-Depth model construction](https://www.linkedin.com/pulse/brief-explanation-3d-density-depth-model-construction-pechnikov/).

## Screenshots

<img src="screenshots/NCubeTopographyBlockSource.jpg" width="400">
<img src="screenshots/NCubeShapefileOnTopographyBlockSource1.jpg" width="400">
<img src="screenshots/NCubeShapefileOnTopographyBlockSource2.jpg" width="400">
<img src="screenshots/NCubeSources.jpg" width="400">

## Used Technologies

We use this set of programming libraries and technologies:

[VTK: The Visualization Toolkit](https://vtk.org/)

[Xarray: N-D labeled arrays and datasets in Python](http://xarray.pydata.org/en/stable/)

[Rasterio: Python access to geospatial raster data](https://rasterio.readthedocs.io/en/latest/)

[GeoPandas: Python geospatial operations on geometric types](http://geopandas.org/)

[Pandas: Python Data Analysis Library](https://pandas.pydata.org/)

[Shapely: Python package for manipulation and analysis of planar geometric objects](https://shapely.readthedocs.io/en/latest/project.html)

[GEOS: Geometry Engine, Open Source](https://trac.osgeo.org/geos/)

[PROJ: generic coordinate transformation software](https://proj.org/)

[Log ASCII Standard (LAS) files in Python](https://lasio.readthedocs.io/en/latest/)

## Authors

Alexey Pechnikov

https://orcid.org/0000-0001-9626-8615 (ORCID)

E-mail: pechnikov@mobigroup.ru
