# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

import sys
import os
sys.path.append(os.path.dirname(__file__))

from paraview.util.vtkAlgorithm import * 

#from NCube import _NCubeGeoDataFrameToTopography, _NCubeGeometryOnTopography
#from NCube import _str, _NCubeGeoDataFrameToTopography, _NCubeGeoDataFrameRowToVTKArrays
from NCube import _NCubeGeometryOnTopography

import xarray as xr
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import box, Point, LineString
from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
import math
import time


#------------------------------------------------------------------------------
# N-Cube Table On Topography Block Source
#
# Create 3D points when X, Y Fields and maybe Z Field defined (overwise Z=0)
# AND
# Create 3D unit vector when Azimuth,Dip Fields defined and Length Field is not defined
#
# OR
#
# Create 3D lines when X,Y,Z,Azimuth,Dip,Length Fields defined
#
#------------------------------------------------------------------------------

@smproxy.source(name="NCubeTableOnTopographyBlockSource",
       label="N-Cube Table On Topography Block Source")
class NCubeTableOnTopographyBlockSource(VTKPythonAlgorithmBase):
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

        self._toponame = None
        self._xcol = None ;# required
        self._ycol = None ;# required
        self._zcol = None ;# optional - topography could be used instead
        self._zunit = 'm'
        self._zinvert = 0
        # depth for earthquakes, etc.
#        self._depthcol = None
#        self._depthunit = 'm'

        self._azcol = None
        self._dipcol = None
        self._dipinvert = 0
        # if it's defined, create segment, overwise ?only create point with vector
        self._lengthcol = None
        self._lengthunit = 'm' ;# required
        self._vectorname = 'vector'


    def RequestData(self, request, inInfo, outInfo):

        if self._tablename is None or self._xcol is None or self._ycol is None:
            return 1

        t0 = time.time()
        #df = pd.read_csv(self._tablename, low_memory=False)

        # load DEM
        dem = dem_extent = dem_crs = None
        if self._toponame is not None:
            #dem = xr.open_rasterio(toponame, chunks=10000000).squeeze()
            dem = xr.open_rasterio(self._toponame).squeeze()
            # dask array can't be processed by this way
            dem.values[dem.values == dem.nodatavals[0]] = np.nan

            # NaN border to easy lookup
            dem.values[0,:]  = np.nan
            dem.values[-1,:] = np.nan
            dem.values[:,0]  = np.nan
            dem.values[:,-1] = np.nan

            #dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())
            #dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None

            #print (dem.values)

#        df = _NCubeDataFrameLoad(self._tablename, self._tablecol,
#            self._xcol, self._ycol, self._zcol, self._zunit,
#            self._tableencoding)
#        if df is None:
#            return

        df = pd.read_csv(self._tablename, encoding=self._tableencoding, low_memory=False)
        print ("len(df)",len(df))
        if (len(df)) == 0:
            return
        
        # TODO
        
        # create point geometry if Length Field is not defined
        xs = df[self._xcol]
        ys = df[self._ycol]
        zs = len(df)*[0]
        if self._zcol is not None:
            zscale = self._scaleunits[self._zunit]
            zs = (1 if self._zinvert == 0 else -1)*zscale*df[self._zcol]

        geom = gpd.GeoSeries(map(Point, zip(xs, ys, zs)))
        crs = {'init' : 'epsg:'+str(self._epsg)} if self._epsg != 0 else None
        print ("crs", crs)
        # crs initialization doesn't work here for Python 2
        df = gpd.GeoDataFrame(df, crs=crs, geometry=geom)
        # for Python 2 only
        df.crs = crs
        print ("df.crs", df.crs)
        # create line segment if Length Field is defined

        # length for unit vector or geometry line
        lengthscale = self._scaleunits[self._lengthunit]
        if self._lengthcol is None:
            length = len(df)*[lengthscale]
        else:
            length = lengthscale*df[self._lengthcol]

        print ("length",length[:5])
        if self._azcol is not None:
            #and self._dipcol is not None:
            # https://github.com/mobigroup/gis-snippets/blob/master/ParaView/ProgrammableFilter/vtkMultiblockDataSet.md
            # https://en.wikipedia.org/wiki/Spherical_coordinate_system
            # Spherical coordinates (r, θ, φ) as often used in mathematics:
            # radial distance r, azimuthal angle θ, and polar angle φ.
            theta = 1./2*math.pi - math.pi*df[self._azcol]/180
            if self._dipcol is not None:
                if self._dipinvert == 0:
                    # 1 for wells (dip<0)
                    phi = math.pi*(90 - df[self._dipcol])/180
                else:
                    # -1 for WSM (dip>0)
                    phi = math.pi*(90 + df[self._dipcol])/180
            else:
                phi = len(df)*[math.pi/2]
            dx = np.round(length*np.sin(phi)*np.cos(theta),10)
            dy = np.round(length*np.sin(phi)*np.sin(theta),10)
            dz = np.round(length*np.cos(phi),10)
            # self._x+row.dx, self._y+row.dy, self._z+row.dz
            df[self._vectorname] = [(0,0,0) if np.any(np.isnan(v)) else v for v in zip(dx,dy,dz)]

            print (df.head()[self._vectorname])
#        # create vector or line geometry
#        if self._lengthcol is None:
#            # create point geometry
#            # create vector
#            if self._azcol is not None and self._dipcol is not None:
#                pass
#        else:
#            # create line geometry

        # group to multiblock
        if self._tablecol is not None:
            df = df.sort_values(self._tablecol).set_index(self._tablecol)
        print (df.head())

        vtk_blocks = _NCubeGeometryOnTopography(df, dem)
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


################
    @smproperty.stringvector(name="Units", information_only="1")
    def GetUnits(self):
        return ["m","km","ft","kft"]

    # possible coordinate columns
    @smproperty.stringvector(name="NumericColumns", information_only="1")
    def GetNumericColumns(self):
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
        if self._tablename is None:
            return []
        # Load file
        df = pd.read_csv(self._tablename, nrows=0)
        cols = sorted(df.columns.values)
        return ['None'] + list(map(str,cols))

################


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

    @smproperty.intvector(name="Table EPSG", default_values=0)
    def SetTableEPSG(self, epsg):
        #epsg = epsg if epsg != 'None' and epsg != 0 else None
        print ("SetTableEPSG", epsg)
        self._epsg = epsg
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



    @smproperty.stringvector(name="Topography File")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["tif", "TIF", "nc"], file_description="GeoTIFF, NetCDF")
    @smdomain.xml(\
        """
        <Hints>
          <PropertyWidgetDecorator type="GenericDecorator"
                                   mode="visibility"
                                   property="TableFileName"
                                   function="boolean_invert" />
        </Hints>
        """)
    def SetTopographyFileName(self, name):
        """Specify filename for the topography file to read."""
        print ("SetTopographyFileName", name)
        name = name if name != 'None' else None
        if self._toponame != name:
            self._toponame = name
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

    @smproperty.stringvector(name="TopElevation", number_of_elements="1")
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

    @smproperty.xml("""
        <IntVectorProperty name="ElevationInvert"
                       command="SetTopZInvert"
                       number_of_elements="1"
                       default_values="0">
        <BooleanDomain name="bool" />
        <Documentation>
            Elevation>0 above the Earth's surface. Use this checkbox to invert (to interpret it as depth).
        </Documentation>
        </IntVectorProperty>
    """)
    def SetTopZInvert(self, value):
        print ("SetTopZInvert", value)
        self._zinvert = value
        #self.Modified()

    @smproperty.stringvector(name="TopElevationUnit", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Units" function="GetUnits"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetTopZUnit(self, unit):
        self._depthunit = unit
        self._zunit = unit
        print("SetTopZUnit", unit)
        self.Modified()



#    @smproperty.stringvector(name="*TopDepth", number_of_elements="1")
#    @smdomain.xml(\
#        """<StringListDomain name="list">
#                <RequiredProperties>
#                    <Property name="NumericColumns" function="GetNumericColumns"/>
#                </RequiredProperties>
#            </StringListDomain>
#        """)
#    def SetTopDepthColumn(self, col):
#        col = col if col != 'None' else None
#        self._depthcol = col
#        print("SetTopDepthColumn", col)
#        self.Modified()
#
#    @smproperty.stringvector(name="*TopDepthUnit", number_of_elements="1")
#    @smdomain.xml(\
#        """<StringListDomain name="list">
#                <RequiredProperties>
#                    <Property name="Units" function="GetUnits"/>
#                </RequiredProperties>
#            </StringListDomain>
#        """)
#    def SetTopDepthUnit(self, unit):
#        self._depthunit = unit
#        self._zunit = unit
#        print("SetTopUnit", unit)
#        self.Modified()




    @smproperty.stringvector(name="DirectionAzimuth", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetDirectionAzimuthhColumn(self, col):
        col = col if col != 'None' else None
        self._azcol = col
        print("SetDirectionAzimuthColumn", col)
        self.Modified()

    @smproperty.stringvector(name="DirectionDip", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetDirectionDipColumn(self, col):
        col = col if col != 'None' else None
        self._dipcol = col
        print("SetDirectionDipColumn", col)
        self.Modified()

    @smproperty.stringvector(name="DirectionLength", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetDirectionLengthColumn(self, col):
        col = col if col != 'None' else None
        self._lengthcol = col
        print("SetDirectionLengthColumn", col)
        self.Modified()

    @smproperty.stringvector(name="DirectionUnit", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Units" function="GetUnits"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetDirectionUnit(self, unit):
        self._lengthunit = unit
        print("SetDirectionUnit", unit)
        self.Modified()

#    @smproperty.intvector(name="Direction Dip Invert", default_values=0, documentation="xxx")
#    @smdomain.xml("""<BooleanDomain name="bool"/>""")
#    def SetDipInvert(self, value):
#        print ("SetDipInvert", value)
#        self._dipinvert = value
#        #self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="DipInvert"
                       command="SetDirectionDipInvert"
                       number_of_elements="1"
                       default_values="0">
        <BooleanDomain name="bool" />
        <Documentation>
            Dip=90° for zenith and dip=-90° for nadir. Use this checkbox to invert.
        </Documentation>
        </IntVectorProperty>
    """)
    def SetDirectionDipInvert(self, value):
        print ("DirectionDipInvert", value)
        self._dipinvert = value
        self.Modified()

    @smproperty.stringvector(name="DirectionVector", default_values="vector")
    def SetDirectionVectorName(self, vname):
        print("SetDirectionVectorName ", vname)
        self._vectorname = vname
        self.Modified()

#    @smproperty.intvector(name="Depth Invert", default_values=0, panel_visibility="advanced")
#    @smdomain.xml(\
#        """<BooleanDomain name="bool"/>
#        <Hints>
#          <PropertyWidgetDecorator type="GenericDecorator"
#                                   mode="visibility"
#                                   property="Depth Unit"
#                                   value="ft" />
#        </Hints>
#        """)
#    def SetYYY(self, epsg):
#        #epsg = epsg if epsg != 'None' and epsg != 0 else None
#        print ("SetYYY", epsg)
#        #self._epsg = epsg
#        #self.Modified()

#    @smproperty.stringvector(name="Output Vector", default_values='vector')
#    def SetVectorName(self, name):
#        print ("SetVectorName", name)
#        self._vectorname = name
#        self.Modified()
