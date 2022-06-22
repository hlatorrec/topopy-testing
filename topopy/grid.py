# -*- coding: utf-8 -*-

# grid.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.1
# February 13, 2018
#
# Last modified 28 September, 2018


from osgeo import gdal
import numpy as np
from scipy import ndimage
from skimage.morphology import reconstruction

# This import statement avoid issues with matplotlib in Mac when using Python not as a Framework
# If matplotlib is not imported, Grid.plot() will not work.
try: 
    import matplotlib.pyplot as plt
    PLT = True
except:
    PLT = False

NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

class PRaster():
    """Defines a Raster object with raster properties and methods. 
    This is a parent class used by other classes to inherit raster properties

    :param path: Path to raster file, creates an empty PRaster when not path is given 
    :type path: str, optional
    """    
    def __init__(self, path=""):        
        if path:
            self._raster = gdal.Open(path)
            raster = self._raster
            if not raster:
                raise FileNotFoundError
            
            self._banda = raster.GetRasterBand(1)
            banda = self._banda
            self._size = (banda.XSize, banda.YSize)
            self._geot = raster.GetGeoTransform()
            self._proj = raster.GetProjection()
        
        else:
            self._size = (1, 1)
            self._geot = (0., 1., 0., 1., 0., -1.)
            self._proj = ""


    def get_size(self):
        """Returns a tuple that contains the size of the grid (XSize, YSize)
        """
        return self._size
    
    def get_dims(self):
        """Returns a tuple that contains the size of the internal array (nrow, ncol)
        """
        return self._size[::-1]
    
    def get_ncells(self):
        """Returns the total number of cells on the Grid
        """
        return self._size[0] * self._size[1]
    
    def get_projection(self):
        """Return a string with the projection of the grid in WKT
        """
        return self._proj
    
    def get_cellsize(self):
        """Returns a tuple with (XCellsize, YCellsize). 
        YCellsize is a negative value
        
        :return: (XCellsize, YCellsize) tuple
        :rtype: tuple
        """
        return self._geot[1], self._geot[5]
    
    def get_geotransform(self):
        """Returns the GeoTransform matrix of the grid. This matrix has the form:
        *(ULx, Cx, Tx, ULy, Ty, Cy)*
        where
        * ULx = Upper-Left X coordinate (upper-left corner of the pixel)
        * ULy = Upper-Left Y coordinate (upper-left corner of the pixel)
        * Cx = X Cellsize
        * Cy = Y Cellsize (negative value)
        * Tx = Rotation on X axis
        * Ty = Rotation on Y axis
        
        :return: GeoTransform matrix
        :rtype: tuple
        """
        return self._geot
    
    def is_inside(self, x, y):
        """Checks if given points are inside the rectangular extent of the raster

        :param x: x coordinate
        :type x: float, list, numpy.ndarray
        :param y: y coordinate
        :type y: float, list, numpy.ndarray
        :return: Boolean or list of booleans (`True` if inside)
        :rtype: boolean, list
        """
        
        row, col = self.xy_2_cell(x, y)
        rowinside = np.logical_and(row >= 0, row < self._size[1])
        colinside = np.logical_and(col >= 0, col < self._size[0])
        inside = np.logical_and(rowinside, colinside)
        return inside
    
    def get_extent(self):
        """Returns a tuple (XMin, XMax, YMin, YMax) with the extension of the Grid
        
        :return: Extension of the grid
        :rtype: tuple
        """
        xmin = self._geot[0]
        xmax = self._geot[0] + self._size[0] * self._geot[1]
        ymin = self._geot[3] + self._size[1] * self._geot[5]
        ymax = self._geot[3]
        
        return (xmin, xmax, ymin, ymax)
    
    def copy_layout(self, grid):
        """Copy all the parameters from another PRaster instance except grid data (and nodata)

        :param pRaster: Raster to copy from 
        :type pRaster: topopy.PRaster  
        """
        self._size = grid.get_size()
        self._geot = grid.get_geotransform()
        self._proj = grid.get_projection()

    def xy_2_cell(self, x, y):
        """Get row and column indexes from XY coordinates
        
        :param x: X coordinate
        :type x: float, list, numpy.ndarray
        :param y: Y coordinates
        :type y: float, list, numpy.ndarray           
        :return: Tuple with (row, col) indices as `numpy.ndarray`
        :rtype: tuple
        """
        x = np.array(x)
        y = np.array(y)       
        row = (self._geot[3] - y) / (self._geot[5] * -1)
        col = (x - self._geot[0]) / self._geot[1]
        return row.astype(np.int32), col.astype(np.int32)

    def cell_2_xy(self, row, col):
        """Get XY coordinates from (row, col) cell indexes
        
        :param row: Row indexes
        :type row: float, list, numpy.ndarray
        :param col: Column indexes
        :type col: float, list, numpy.ndarray
        :return: Tuple with (x, y) coordinates as `np.ndarray`
        :rtype: tuple
        """       
        row = np.array(row)
        col = np.array(col)
        x = self._geot[0] + self._geot[1] * col + self._geot[1] / 2
        y = self._geot[3] + self._geot[5] * row + self._geot[5] / 2
        return x, y
    
    def ind_2_cell(self, ind):
        """
        Get row col indexes from cells linear indexes (row-major, C-style)
        
        :param ind: Linear indexes
        :type ind: float, list, np.ndarray
        :return: Tuple with (row, col) indices as `numpy.ndarray`
        :rtype: tuple
        """
        return np.unravel_index(ind, self.get_dims()) 
    
    def cell_2_ind(self, row, col):
        """Get cell linear indexes from row and column indexes
        
        :param row: Row indexes
        :type row: float, list, numpy.ndarray
        :param col: Column indexes
        :type col: float, list, numpy.ndarray
        :return: Array containing linear indices (row-major, C-style)
        :rtype: np.ndarray
        """
        return np.ravel_multi_index((row, col), self.get_dims())
    
class Grid(PRaster):
    """Class to manipulate rasters
    
    :param path: Path to raster file    
    :type path: str, optional
    :param band: Raster band to be opened (usually doesn't need to be modified)
    :type band: int, optional
    """
    def __init__(self, path="", band=1):    
        # Elements inherited from PRaster.__init__
        super().__init__(path)

        # New elements of Grid        
        if path: 
            banda = self._banda
            self._nodata = banda.GetNoDataValue()
            self._array = banda.ReadAsArray()
            self._tipo = str(self._array.dtype)
                  
        else:
            self._nodata = None
            self._array = np.array([[0]], dtype=np.float)
            self._tipo = str(self._array.dtype)
           
    def set_array(self, array):
        """Set the data array for the current Grid object. 
        If the current Grid is an empty Grid [`get_size()` = (1, 1)], any input array is valid
        If The current Grid is not an empty Grid, the input array should match Grid dimensions
        
        :param array: Numpy array containing the data
        :type array: numpy.ndarray
        """
        # If the Grid is an empty Grid, any array is valid
        if self._size == (1, 1):       
            self._size = (array.shape[::-1])    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
        # If the Grid is not an empty Grid, input array shape must coincide with internal array
        elif array.shape == self.get_dims():    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
        else:
            return 0
       
    def read_array(self, ascopy=False):
        """Reads the internal array of the Grid instace
        
        :param ascopy: If `True`, the returned array is a memory view of the original array
        :type ascopy: bool, optional
        :return: Internal array of the current Grid object
        :rtype: numpy.ndarray
        """
        if ascopy:
            return np.copy(self._array)
        else:
            return self._array
    
    def find(self):
        """Find the non-zero elements in the array.
        
        :return: Tuple of arrays with row and col positions
        :rtype: tuple
        """
        return np.where(self._array > 0)
    
    def max(self):
        """Get maximum value of the Grid
        
        :return: Max value
        :rtype: float
        """
        datapos = np.where(self._array != self._nodata)
        return np.max(self._array[datapos])
    
    def min(self):
        """Get minimun value of the Grid
        
        :return: Min value
        :rtype: float
        """
        datapos = np.where(self._array != self._nodata)
        return np.min(self._array[datapos])
    
    def mean(self):
        """Get mean value of the Grid
        
        :return: Mean value
        :rtype: float
        """
        datapos = np.where(self._array != self._nodata)
        return np.mean(self._array[datapos])
    
    def set_value(self, row, col, value):
        """Sets value for the specified cells at (row, col)
        
        :param row: Row index
        :type row: int, numpy.ndarray
        :param col: Column index
        :type col: int, numpy.ndarray
        :param value: Value to set on the cell
        :type value: float
        """
        self._array[row, col] = value
    
    def get_value(self, row, col):
        """Gets value of the specified cells at (row, col)
        
        :param row: Row index
        :type row: int, numpy.ndarray
        :param col: Column index
        :type col: int, numpy.ndarray
        :return: Cell values
        :rtype: float, numpy.ndarray
        """
        return self._array[row, col]
    
    def get_nodata(self):
        """Get value for NoData cells in the grid. This value could be `None` if NoData
        is not defined in the grid.
        
        :return: NoData value
        :rtype: float
        """
        return self._nodata
    
    def set_nodata(self, value):
        """Sets NoData value for the Grid
        
        :param value: NoData value
        :type value: float
        """
        # If nodata wasn't stabished, we set up the new value
        if self._nodata is None:
            self._nodata = value
        # If value == None, we are removing nodata value
        elif value is None:
            self._nodata = None
        # If grid had a nodata defined, we set up the new value
        # And change old nodata values with the new value
        else:
            self._array[self.get_nodata_pos()] = value
            self._nodata = value
            
    def get_nodata_pos(self):
        """Get positions of NoData values as a tuple of two arrays (rows, columns)
        
        :return: Indices locating NoData values of the Grid
        :rtype: tuple
        """
        if self._nodata is None:
            return (np.array([], np.int), np.array([], np.int))
        else:
            return np.where(self._array == self._nodata)
        
    def nan_2_nodata(self):
        """Changes `NaN` values to NoData (if Grid NoData is defined).
        """
        if self._nodata is None:
            return
        idx = np.isnan(self._array)
        self._array[idx] = self._nodata
    
    def values_2_nodata(self, value):
        """Change specific values to NoData (if Grid NoData is defined).
        
        :param value: Values to change
        :type value: float
        """
        if self._nodata is None:
            return
        
        if type(value) == int or type(value)==float:
            ind = np.where(self._array==value)
            self._array[ind] = self._nodata
        else:
            for val in value:
                ind = np.where(self._array == val)
                self._array[ind] = self._nodata        
    
    def plot(self, ax=None):
        """Plots the grid in a new Axes or an existing one
        
        :param ax: Defaults to `plt.imshow` if not defined
        :type ax: matplotlib.Axe, optional
        """
        if not PLT:
            return
        
        arr = self.read_array()
        if self._nodata:
            mask = self._array == self._nodata
            arr = np.ma.array (self._array, mask = mask)
        
        if ax:
            ax.imshow(arr)
        else:
            plt.imshow(arr)
    
    def is_inside(self, x, y, NoData=True):
        """Checks if points are inside the raster
        
        :param x: x coordinate
        :type x: float, list, numpy.ndarray
        :param y: y coordinate
        :type y: float, list, numpy.ndarray
        :param NoData: Flag to consider NoData values, default `True` (NoData values are outside of the raster)
        :type NoData: bool, optional
        :return: Boolean or list of booleans, inside is `True` and outside is `False`
        :rtype: bool, list
        """
        inside = super().is_inside(x, y)

        if NoData:
            pos = np.where(inside)
            x = x[inside]
            y = y[inside]
            row, col = self.xy_2_cell(x, y)
            in_nodata = self.get_value(row, col) != self.get_nodata()
            inside[pos] = in_nodata
        
        return inside
    
    def copy(self):    
        """Copy the Grid
        
        :return: Copy of the Grid
        :rtype: topopy.Grid
        """
        newgrid = Grid()
        newgrid.copy_layout(self)
        newgrid._array = np.copy(self._array)
        newgrid._tipo = self._tipo
        newgrid._nodata = self._nodata
        return newgrid
    
    def save(self, path):
        """Saves the Grid onto disc
        
        :param path: Path to save the Grid at
        :type path: str
        """
        # Check if the type of the internal array is compatible with gdal
        if str(self._tipo) not in NTYPES.keys():
            return 0
        else:
            tipo = NTYPES[str(self._tipo)]
        
        # Prepare driver to write raster
        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._size[0], self._size[1], 1, tipo)
        if not raster:
            return
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)
        
        if not self._nodata is None:
            raster.GetRasterBand(1).SetNoDataValue(self._nodata)
        raster.GetRasterBand(1).WriteArray(self._array)

class DEM(Grid):
    
    def __init__(self, path="", band=1):
        
        # Call to Grid.__init__ method
        super().__init__(path, band)
        # Change NoData values to -9999.
        self.set_nodata(-9999.)
        
    def copy(self):
        """Copy DEM
        
        :return: DEM object
        :rtype: topopy.DEM
        """
        newgrid = DEM()
        newgrid.copy_layout(self)
        newgrid._array = np.copy(self._array)
        newgrid._tipo = self._tipo
        newgrid._nodata = self._nodata
        return newgrid
        
    def identify_flats(self, nodata=True, as_array=False):
        """This functions returns two `topopy.Grid` (or `numpy.ndarray`) with flats and sills. 
        Flats are defined as cells without downward neighboring cells. Sills are cells where 
        flat regions spill over into lower terrain. It the DEM has NoData values, those will be maintained
        in output grids.
        
        References:
        -----------
        This algoritm is adapted from identifyflats.m by Wolfgang Schwanghart 
        (version of 17. August, 2017) included in TopoToolbox matlab codes.
        
        Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
        for topographic analysis. Environ. Model. Softw. 25, 770–781. 
        https://doi.org/10.1016/j.envsoft.2009.12.002
        
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        
        :param nodata: Flag to keep or replace NoData values (default `True`).
            If `True`, values are kept. If `False`, values are repolaced by zeros
        :type nodata: bool, optional
        :param as_array: Flag to chose between `numpy.ndarray` and `topopy.Grid` outputs (default `False`).
            If `True`, outputs are `numpy.array`. If `False`, outputs are `topopy.Grid`
        :type as_array: bool, optional
        :return: Tuple with output Grids (flats, sills)
        :rtype: tuple
        """

        z_arr = np.copy(self._array)
        
        # Change nodata to very low values
        nodata_ids = self.get_nodata_pos()
        z_arr[nodata_ids] = -9999
        
        footprint = np.ones((3, 3), dtype=np.int8)
        # Identify flats throught a image binary erosion
        # Flats will be True where cells don't have lower neighbors
        flats = ndimage.morphology.grey_erosion(z_arr, footprint=footprint) == z_arr
        
        # Remove flats from the borders
        flats[0,:] = flats[-1,:] = flats[:,0] = flats[:,-1] = False
        
        # Remove flats for nodata values and cells bordering them
        flats[nodata_ids] = False
        auxmat = np.zeros(flats.shape, dtype="bool")
        auxmat[nodata_ids] = True
        nodata_bord = ndimage.morphology.grey_dilation(auxmat, footprint=footprint)
        flats[nodata_bord] = False
        
        # Identify sills
        sills = np.empty(z_arr.shape)
        sills.fill(-9999.)
        sills[flats] = z_arr[flats]
        aux_dil = ndimage.morphology.grey_dilation(sills, footprint=footprint)
        sills = np.logical_and(aux_dil == z_arr, np.logical_not(flats))
        sills[nodata_ids] = False
        
        # Prepare outputs
        if as_array:
            return flats, sills
        else:
            res = []
            for arr in [flats, sills]:
                grid = Grid()
                grid.copy_layout(self)
                grid.set_nodata(-1)
                grid._tipo = 'int8'
                grid.set_array(arr.astype(np.int8))
                
                if nodata:
                    grid._array[nodata_ids] =  -1
                else:
                    grid._array[nodata_ids] =  0
                    grid.set_nodata(None)
            
                res.append(grid)
    
            return res

    def fill_sinks(self, as_array=False):
        """Fill sinks in a DEM using scikit-image reconstruction algorithm
        
        :param as_array: Flag to return DEM as `numpy.array` (default `False`) instead of `topopy.DEM`
        :type as_array: bool, optional
        :return: Filled DEM
        :rtype: topopy.DEM
        """
        # Get the seed to start the fill process
        seed = np.copy(self._array)
        seed[1:-1, 1:-1] = self._array.max()

        # Fill the DEM        
        nodata_pos = self.get_nodata_pos()
        filled = reconstruction(seed, self._array, 'erosion')
        filled = filled.astype(self._array.dtype)
        filled[nodata_pos] = self._nodata
        
        if as_array:
            # Return filled DEM as numpy.ndarray
            return filled
        else:
            # Return filled DEM as topopy.DEM
            filled_dem = self.copy()
            filled_dem.set_array(filled)
            return filled_dem
    
    def fill_sinks2(self, four_way=False):
        """Fill sinks method adapted from fill depressions/sinks in floating point array
        
        This algorithm has been adapted (with minor modifications) from the 
        Charles Morton slow fill algorithm (with ndimage and python 3 it was not slow
        at all). 
        
        References
        ----------
        Soile, P., Vogt, J., and Colombo, R., 2003. Carving and Adaptive
        Drainage Enforcement of Grid Digital Elevation Models.
        Water Resources Research, 39(12), 1366
        
        Soille, P., 1999. Morphological Image Analysis: Principles and
        Applications, Springer-Verlag, pp. 173-174
        
        :param four_way: Chose between 4 (`True`) adjacent cells or 8 (`False`) adjacent cells
        :type four_way: bool, optional
        :return: Filled DEM
        :rtype: topopy.DEM 
        """
        # Change nan values to a very low value
        copyarr = np.copy(self._array)
        nodata_pos = self.get_nodata_pos()
        copyarr[nodata_pos] = -9999.
        
        # Set h_max to a value larger than the array maximum to ensure
        #   that the while loop will terminate
        h_max = copyarr.max() + 100
    
        # Build mask of cells with data not on the edge of the image
        # Use 3x3 square Structuring element
        inside_mask = ndimage.morphology.binary_erosion(np.isfinite(copyarr),np.ones((3,3), bool))
    
        # Initialize output array as max value test_array except edges
        output_array = np.copy(copyarr)
        output_array[inside_mask] = h_max
    
        # Array for storing previous iteration
        output_old_array = np.copy(copyarr)
        output_old_array[:] = 0
    
        # Cross structuring element
        if four_way:
            el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
        else:
            el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]], bool)
    
        # Iterate until marker array doesn't change
        while not np.array_equal(output_old_array, output_array):
            output_old_array = np.copy(output_array)
            output_array = np.maximum(
                copyarr,
                ndimage.grey_erosion(output_array, size=(3, 3), footprint=el))

        # Put back nodata values and change type
        if self._nodata:
            output_array[nodata_pos] = self._nodata
        # Create output filled DEM
        filled_dem = DEM()
        filled_dem.copy_layout(self)
        filled_dem.set_array(output_array)
        filled_dem.set_nodata(self._nodata)
        return filled_dem

class Basin(DEM):
    """Class to manipulate drainage basins. The object is basically a DEM with NoData
    values in cells outside of the drainage basin.
    
    :param dem: DEM instance. If dem is `str` and basin is `None`, it will load from the specified path
    :type dem: str, topopy.DEM
    :param basin: Drainage basin. If `None`, DEM is loaded as a basin. Needs same dimensions and cellsize as the input DEM
    :type basin: None, str, topopy.Grid, optional
    :param idx: Value of the basin cells
    :type idx: int, optional       
    """
    def __init__(self, dem, basin=None, idx=1):
        # If basin is None, the DEM is already a basin and we load it
        if basin is None:
            # Call to DEM.__init__ method
            super().__init__(dem)
            # Change NoData values to -9999.
            self.set_nodata(-9999.)
            return
        elif type(basin) is str:
            basingrid = Grid(basin)           
        elif type(basin) is Grid:
            basingrid = basin
        
        #dem = DEM(dem, 1)
        basin = np.where(basingrid.read_array()==idx, 1, 0)
        
        if basingrid.get_size() != dem.get_size():
            raise GridError("ERROR. DEM and basin grids have different dimensions!")

        # Get limits for the input basin
        c1 = basin.max(axis=0).argmax()
        r1 = basin.max(axis=1).argmax()
        c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax()
        r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax()
        
        # Cut basin and dem by those limits
        basin_cl = basin[r1:r2, c1:c2]
        dem_cl = dem.read_array()[r1:r2, c1:c2]

        # Create Grid
        self._size = (basin_cl.shape[1], basin_cl.shape[0])
        self._dims = (basin_cl.shape[0], basin_cl.shape[1])
        geot = dem._geot
        ULx = geot[0] + geot[1] * c1
        ULy = geot[3] + geot[5] * r1
        self._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
        self._cellsize = (geot[1], geot[5])
        self._proj = dem._proj
        self._ncells = basin_cl.size
        self._nodata = dem._nodata
        self._tipo = dem._tipo
        arr = np.where(basin_cl > 0, dem_cl, dem._nodata)
        self._array = arr.astype(dem._tipo)
        
    def copy(self):    
        """Get a copy of the Basin
        
        :return: Basin object
        :rtype: topopy.Basin
        """
        newgrid = Basin()
        newgrid.copy_layout(self)
        newgrid._array = np.copy(self._array)
        newgrid._tipo = self._tipo
        newgrid._nodata = self._nodata
        return newgrid

class GridError(Exception):
    pass