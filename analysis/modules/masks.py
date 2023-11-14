import numpy as np
from gammapy.maps import Map

# Apply mask to map
def apply_mask(map, mask_map):
    if not isinstance(mask_map, Map) or not isinstance(map, Map):
        raise ValueError("Both input maps must be gammapy.maps.Map objects.")

    map.data = map.data - (map.data * mask_map)
    # Compute the median of the non-zero elements
    m = np.median(map.data[map.data > 0])
    # Assign the median to the zero elements 
    map.data[map.data == 0] = m
    
    return map