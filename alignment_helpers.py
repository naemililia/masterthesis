import numpy as np
from pyproj import CRS, Transformer
import re
import geopandas as gpd
from shapely.geometry import Point
import utm
import math 

def get_coordinates_of_corners(transform_sen, array_sen: np.ndarray) -> dict:
    """
    
    This function returns the coordinates of the 4 corners of the unlabeled sentinel image.
    
    Careful: this only works when you read only 1 band from the sentinel image

    TODO: think about whether it makes sense, to include the UTM to WGS projection in here already, or if i want to keep it explicit

    Args:
        transform_sen (_type_): affine transformation matrix of the unlabeled sentinel image
        array_sen (np.ndarray): array of the unlabeled sentinel image 

    Returns:
        dict: the coordinates of the 4 corners in format (lon, lat)
    """
    
    max_x,max_y = array_sen.shape
    
    top_left_lon, top_left_lat = transform_sen * (0,0)
    bottom_left_lon, bottom_left_lat = transform_sen * (0, max_y-1)
    bottom_right_lon, bottom_right_lat = transform_sen * (max_x-1, max_y-1)
    top_right_lon, top_right_lat = transform_sen * (max_x-1, 0)
    
    corner_coordinates = {
    'top_left' : (top_left_lon, top_left_lat),
    'top_right' : (top_right_lon, top_right_lat),
    'bottom_left' : (bottom_left_lon, bottom_left_lat),
    'bottom_right' : (bottom_right_lon, bottom_right_lat)
    }
    
    return corner_coordinates


def get_x_and_y_pixel(target_lat: float, target_lon: float, transform_aux):
    """
    Description:
    Given the coordinates of a point in an unlabeled satellite image, it finds the corresponding pixels of the same location in the auxiliary datasource.
    
    Example functionality:
    (target_lon - transform_aux[2]) = difference of target lon and the lon that is at the top left corner --> degrees of longtitude that need to be 'crossed' to reach the taregt lon
    ! minus is correct here, because: the target lon will always be equal or further right than the top left corner, making the value always positive, no matter whether the location is on the east (negative lon) or west (positive lon) or on the edge of east and west
    
    this difference is divided by transform_aux[0], which is a value describing how many degrees are 'crossed' per pixel
    
    the value x_pixel is therefore the number of pixels needed to move to the right of the left side of the image to cross the amount of degrees that the target lon is away from the lon at the left side
    
    Because we are working with whole pixels, this number is rounded. Because of that rounding, we are introducing a slight margin of error
    
    Inputs:
    target_lat: latitude of the location i want to find the corresponding pixel in the auxiliary info for
    target_lon: longitude of the location i want to find the correspnding pixel in the auxiliary info for
    transform_aux: affine trasnformation matrix of the auxiliary information (this needs to be in same CRS as the target_lon & target_lat)
    """
    x_pixel = round((target_lon - transform_aux[2]) / transform_aux[0])
    y_pixel = round((target_lat - transform_aux[5]) / transform_aux[4])
    
    return x_pixel, y_pixel

def get_corner_pixels(corner_coordinates, transform_aux):
    """
    This function is used to spatially align two satellite images. 
    Loops through the dictionary storing the corner coordinates of the satelite image, and calls get_x_and_y_pixel on each

    Args:
        corner_coordinates (dict): WGS84 coordinates of corners of unlabeled staelite image
        transform_aux: affine transformation of auxiliary data source

    Returns:
        dict: 4 keys, 1 for each corner, each carrying (x,y) pixel coordinates
    """
    corner_pixels = {}

    for key, value in corner_coordinates.items():
        x,y = get_x_and_y_pixel(value[1], value[0], transform_aux)
        corner_pixels[key] = (x,y)
    
    return corner_pixels


def get_label_from_pixel(corner_pixels, array_aux):
    """
    slice the aux array according to the pixel values which enircle the entire area of interest, i.e. the unlabeled image --> create the new weak label

    Args:
        corner_pixels (dict): carries the pixel coorinates in teh aux info corresponding to the geocoordinates in teh unlabeled satelite image 
        array_aux (array): array belonging to the auxiliary information (has to be only 1 layer/band (cant be RGB))

    Returns:
        array: array sliced at the pixel coordinates
    """
    return array_aux[corner_pixels['top_left'][1]: corner_pixels['bottom_left'][1], corner_pixels['top_left'][0]:corner_pixels['top_right'][0]]


def find_UTM_zone_eurosat(crs_string):
    """
    Identifies utm zone (zone nr and hemisphere, ex: '32N') from teh crs

    Args:
        crs_string (str): the entire coordinate reference system saved as a string (done during opening of geotif file)

    Returns:
        zone_nr (str): the UTM zone
        hemisphere (str): 'N' or 'S', the UTM hemisphere
    """
    # describe the pattern in which the UTM zone is described in teh eurosat images
    pattern = r'UTM\szone\s\d{1,2}[NS]'

    # Search for the pattern in the crs_string
    match = re.search(pattern, crs_string)

    # identify the components
    zone = match.group().split()[2]
    zone_nr = zone[:2]
    hemisphere = zone[2:]
    
    return zone_nr, hemisphere

def get_coordinates_from_utm(x_coord, y_coord, zone_nr, hemisphere):
    
    
    lat, lon = utm.to_latlon(x_coord, y_coord, zone_nr, hemisphere)
    
    # Option with pyproj commented out below
    
    # # Define CRS for UTM zone and WGS84
    # utm_crs = CRS(proj='utm', zone=zone_nr, ellps='WGS84', datum='WGS84')
    # wgs84_crs = CRS(proj='latlong', datum='WGS84')

    # # Create a transformer
    # transformer = Transformer.from_crs(utm_crs, wgs84_crs)

    # # Example UTM coordinates (easting, northing) in zone 
    # utm_easting = x_coord
    # utm_northing = y_coord

    # # Convert UTM to WGS84
    # lon, lat = transformer.transform(utm_easting, utm_northing)
    
    return (lon, lat)

def get_coordinates_of_corners_from_UTM(zone_nr, corner_coordinates_euro, hemisphere):

    wgs_coordinates_euro = {}

    for key, value in corner_coordinates_euro.items():
        lon, lat = get_coordinates_from_utm(value[0], value[1], zone_nr, hemisphere)
        wgs_coordinates_euro [key] = (lon,lat)
        
    return wgs_coordinates_euro

def get_BUILT_polygon_encircling_unlabeled_staelite_image(shapefile_of_tiles, corner_coordinates_dict: dict):
    """
    
    this function will return the tile_id of the BUILT tiles which encircles the entire unlabeled satelite image 
    
    if the image does relate to a single tile, return the tile_id
    
    if the image is spread out over a multitude of tiles, for now, return an error / skip that one
    TODO: handle the case where it doesnt relate to a single on

    Args:
        corner_coordinates_dict (dict): dictionary that results form the function get_coordinates_of_corners
    """
    
    tile_ids = []
    for key, value in corner_coordinates_dict.items():
        search_point = Point(value[0], value[1])
        # Find the row of shapefile where the polygon encloses the specific coordinate
        # This creates a boolean series where True represents polygons containing the point
        contains_point_series = shapefile_of_tiles['geometry'].contains(search_point)
        # Filter the shapefile to only include rows where the polygon contains the point
        enclosing_polygons = shapefile_of_tiles[contains_point_series]
        
        # I just want the first row that matches (assuming at least one match exists), as the polygons do not overlap
        
        # TODO: figure out how to check if the length of that enclosing polygon is larger than 1, and then ony
        if len(enclosing_polygons) == 1:
            tile_id_first_enclosing_polygon = enclosing_polygons.iloc[0]['tile_id']
            tile_ids.append(tile_id_first_enclosing_polygon)
        elif len(enclosing_polygons) == 0:
            print('no tile enclosing the satelite image found')
            
    # now i need to check if all the returned tiles are the same
    unique_entries = set(tile_ids)
    number_of_unique_entries = len(unique_entries)
    if number_of_unique_entries == 1:
        print(f'{tile_ids[0]} is the only tile to download')
        return tile_ids[0]
    else:
        print('there is more than one tile needed, the image is spread out over amultitude of tiles')

def next_biggest_multiple_of_20(number):
    """
    Returns the next biggest multiple of 20 of the absolute value of the provided number.
    """
    if number % 20 == 0:
        return number + 20
    else:
        return number + (20 - number % 20)
    
def next_smallest_multiple_of_20(number):
    """
    Returns the next biggest multiple of 20 of the absolute value of the provided number.
    """
    if number % 20 == 0:
        return number
    else:
        return number - (number % 20)


def find_correct_land_cover_tile(corner_coordinates_sat):
    
    """
    this function identifies the tile id of the land cover auxiliary data source in which the satelite image is located.
    
    does this by taking the top left coordinates of the satelite image and finding the nearest 20 degree grid mark

    Returns:
        string: the string with which the tile files start. Entire file name then an be created like so: 
                f'{tile_id_lc}_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif'
    """
    nbm_lon = next_smallest_multiple_of_20(corner_coordinates_sat['top_left'][0])
    nsm_lat = next_biggest_multiple_of_20(corner_coordinates_sat['top_left'][1])
    
    if corner_coordinates_sat['top_left'][0] <= 0:
        lon = 'W'
    else: lon = 'E'
    
    if corner_coordinates_sat['top_left'][1] <= 0:
        lat = 'S'
    else: lat = 'N'
    
    if nbm_lon == 0:
        lon = 'E'
        
    if nsm_lat == 0:
        lat = 'N'
        
    longt = lon+str(abs(int(nbm_lon)))
    latdt = lat+str(abs(int(nsm_lat)))

    # wenn das jetzt nur 3 stellen hat, dann inserte eine 0 an stelle 2

    if len(longt) == 3:
        longt = longt[:1] + '0' + longt[1:]
        
    elif len(longt) == 2:
        longt = longt[:1] + '00' + longt[1:]
        
    if len(latdt) == 2:
        latdt = latdt[:1] + '0' + latdt[1:]
        
    tile_id_lc = longt+latdt
    
        
    return tile_id_lc

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    
    Basically find the distance between two coordinates.
    Used to check how misaligned the label and the aux info is
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r