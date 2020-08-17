import requests
import random

def get_altitdue(lat, lon):
	"""
	returning the altitude of the point, uses free mapquest data that is limited in the amount of calls per month, it uses Nir BAMBA Benmoshe key

	param lat - the path where we save the stl
	param lon - the width of the domain in the x direction

	return:
	altitude - meters above sea level
	"""

	resp = requests.get('http://open.mapquestapi.com/elevation/v1/profile?key=D5z9RSebQJLbUs4bohANIB4TzJdbvyvm&shapeFormat=raw&latLngCollection='+str(lat)+','+str(lon))


	height = resp.json()['elevationProfile'][0]['height']

	return height

if __name__ == "__main__":
    lon = random.randint(35750, 35800) / 1000.0 # Hermon
    lat = random.randint(33250, 33800) / 1000.0
    alt = get_altitdue(lat,lon)
    print("the altitude at position: ",lat,lon," is ", alt)


