# -*- coding: utf-8 -*-
# Copyright (c) 2020 Alexey Pechnikov. All rights reserved.
# https://orcid.org/0000-0001-9626-8615 (ORCID)
# pechnikov@mobigroup.ru (email)
# License: http://opensource.org/licenses/MIT

from paraview.util.vtkAlgorithm import * 

from NCube import _NCubeGeoDataFrameToTopography

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
        self._tablename = None ;# required
        self._tablecol = None
        self._tableencoding = None
        self._toponame = None
        self._xcol = None ;# required
        self._ycol = None ;# required
        self._zcol = None
        self._zunit = 'm' ;# required
        self._epsg = None

        self._azcol = None
        self._dipcol = None
        # if it's defined, create segment, overwise ?only create point with vector
        self._lengthcol = None
        self._lengthunit = 'm' ;# required
        self._vectorname = 'vector'


    def RequestData(self, request, inInfo, outInfo):
        import geopandas as gpd
        import pandas as pd
        from shapely.geometry import Point, LineString
        from vtk import vtkPolyData, vtkAppendPolyData, vtkCompositeDataSet, vtkMultiBlockDataSet
        import time

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

            dem_extent = box(dem.x.min(),dem.y.min(),dem.x.max(),dem.y.max())
            dem_crs = dem.crs if 'crs' in dem.attrs.keys() else None

            #print (dem.values)

#        df = _NCubeDataFrameLoad(self._tablename, self._tablecol,
#            self._xcol, self._ycol, self._zcol, self._zunit,
#            self._tableencoding)
#        if df is None:
#            return

        df = pd.read_csv(self._tablename, encoding=self._tableencoding, low_memory=False)
        if (len(df)) == 0:
            return
        if self._tablecol is not None:
            df = df.sort_values(self._tablecol).set_index(self._tablecol)
        #print ("shapecol",shapecol)
        
        # TODO
        
        # create point geometry if Length Field is not defined
        geom = gpd.GeoSeries(map(Point, zip(df[self._xcol], df[self._ycol], df[self._zcol] if self._zcol is not None else len(df)*[0])))
        #df.drop(['x', 'y'], axis=1, inplace=True)
        df = gpd.GeoDataFrame(df, geometry=geom)
        print (df.head())
        # create line segment if Length Field is defined
        
        if dem is not None:
            df = _NCubeGeoDataFrameToTopography(df, dem_extent, dem_crs)

        groups = df.index.unique() ;#[11454:11455]

#        vtk_blocks = _NCubeGeometryOnTopography(self._shapename, self._toponame, self._shapecol, self._shapeencoding)
#        if vtk_blocks is None or vtk_blocks == []:
#            t1 = time.time()
#            print ("t1-t0", t1-t0)
#            return 1
#        print ("vtk_blocks", len(vtk_blocks))
#        mb = vtkMultiBlockDataSet.GetData(outInfo, 0)
#        mb.SetNumberOfBlocks(len(vtk_blocks))
#        rowidx = 0
#        for (label, polyData) in vtk_blocks:
#            #print (rowidx, label)
#            mb.SetBlock( rowidx, polyData )
#            mb.GetMetaData( rowidx ).Set( vtkCompositeDataSet.NAME(), label)
#            rowidx += 1
        t1 = time.time()
        print ("t1-t0", t1-t0)

        return 1

    @smproperty.stringvector(name="Table File Name")
    @smdomain.filelist()
    @smhint.filechooser(extensions=["csv"], file_description="CSV")
    def SetShapeFileName(self, name):
        """Specify filename for the table to read."""
        print ("SetCSVFileName", name)
        name = name if name != 'None' else None
        if self._tablename != name:
            self._tablename = name
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

    # Multiblock labels
    @smproperty.stringvector(name="Columns", information_only="1")
    def GetColumns(self):
        if self._tablename is None:
            return []
        # Load file
        import pandas as pd
        df = pd.read_csv(self._tablename, nrows=0)
        cols = sorted(df.columns.values)
        return ['None'] + list(map(str,cols))
    @smproperty.stringvector(name="Group by (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Columns" function="GetColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetGroup(self, col):
        col = col if col != 'None' else None
        self._tablecol = col
        print("SetGroup", col)
        self.Modified()

    # possible coordinate columns
    @smproperty.stringvector(name="NumericColumns", information_only="1")
    def GetNumericColumns(self):
        print ("GetColumns 0")
        if self._tablename is None:
            return []
        # Load file
        import pandas as pd
        df = pd.read_csv(self._tablename, nrows=1)
        cols = sorted(df._get_numeric_data().columns)
        print ("GetColumns", cols)
        return ['None'] + list(map(str,cols))

    @smproperty.stringvector(name="Easting", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetXColumn(self, col):
        col = col if col != 'None' else None
        self._xcol = col
        print("SetXColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Northing", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetYColumn(self, col):
        col = col if col != 'None' else None
        self._ycol = col
        print("SetYColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Elevation (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetZColumn(self, col):
        col = col if col != 'None' else None
        self._zcol = col
        print("SetZColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Azimuth (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetAzimuthColumn(self, col):
        col = col if col != 'None' else None
        self._azcol = col
        print("SetAzColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Dip (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetDipColumn(self, col):
        col = col if col != 'None' else None
        self._dipcol = col
        print("SetDipColumn", col)
        self.Modified()

    @smproperty.stringvector(name="Length (optional)", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="NumericColumns" function="GetNumericColumns"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetLengthColumn(self, col):
        col = col if col != 'None' else None
        self._lengthcol = col
        print("SetLengthColumn", col)
        self.Modified()


    @smproperty.stringvector(name="Units", information_only="1")
    def GetUnits(self):
        return ["m","km","ft","kft"]

    @smproperty.stringvector(name="Length Unit", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Units" function="GetUnits"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetLengthUnit(self, unit):
        self._lengthunit = unit
        print("SetLengthUnit", unit)
        self.Modified()

    @smproperty.stringvector(name="Z Unit", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="Units" function="GetUnits"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetZUnit(self, unit):
        self._zunit = unit
        print("SetZUnit", unit)
        self.Modified()

    @smproperty.intvector(name=" EPSG (optional)", default_values=0)
    def SetEPSG(self, epsg):
        self._epsg = epsg
        self.Modified()

    @smproperty.stringvector(name="Output Vector", default_values='vector')
    def SetVectorName(self, name):
        print ("SetVectorName", name)
        self._vectorname = name
        self.Modified()
