from .image import datalayer as ImageDatalayer
from .topography import datalayer as TopographyDatalayer
from .buildings import datalayer as BuildingsDatalayer
topography = TopographyDatalayer("Topography")
buildings = BuildingsDatalayer("Buildings")
image = ImageDatalayer('imageLocations')