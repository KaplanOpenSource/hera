from .demography import datalayer as demographyDatalayer
demography = demographyDatalayer("Demography")

from .locations.image import datalayer as ImageDatalayer
from .locations.topography import datalayer as TopographyDatalayer
from .locations.buildings import datalayer as BuildingsDatalayer
topography = TopographyDatalayer("Topography")
buildings = BuildingsDatalayer("Buildings")
#image = ImageDatalayer('imageLocations')