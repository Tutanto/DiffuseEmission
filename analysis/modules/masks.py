import numpy as np
from gammapy.maps import Map

# Apply mask to counts
def mask_counts(counts, mask_map):
    if not isinstance(mask_map, Map) or not isinstance(counts, Map):
        raise ValueError("Both input maps must be gammapy.maps.Map objects.")

    counts.data = counts.data - (counts.data * mask_map)
    # Compute the median of the non-zero elements
    m = np.median(counts.data[counts.data > 0])
    # Assign the median to the zero elements 
    counts.data[counts.data == 0] = m
    
    return counts

# Apply mask to excess
def mask_excess(excess, mask_map):
    if not isinstance(mask_map, Map) or not isinstance(excess, Map):
        raise ValueError("Both input maps must be gammapy.maps.Map objects.")
    excess.data = excess.data - (excess.data * mask_map)
    # Compute the median of the non-zero elements
    m = np.median(excess.data[excess.data > 0])
    # Assign the median to the zero elements 
    excess.data[excess.data == 0] = m 
    
    return excess