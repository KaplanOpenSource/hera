import numpy
import math
import pandas


class spatialInterpolate():

    def interp(self, point, stations, topography=None, dx=20, dy=20, C=1000, D=5, Hsl=100, b=150):
        """
        Interpolate the values of a variable by its values in different stations.
        params:
            point: The point ([x,y,z])
            topography: The height of the topography at the point
            stations: list of information from stations. Each value in the list
                      is of the structure [latitude, longitude, elevation,[values of variable]]
            dx,dy: Measures of each cell
            C: The value of parameter C, which represent the measure of influence between layers in the atmospheres high above the ground.
            D: The value of parameter D, which represent the measure of influence between layers in the atmospheres near the ground.
            Hsl: The value of parameter Hsl, which effect the height above ground in which the behavor changes from D to C.
            b: The value of parameter b, which effect the slope of the change from D to C.
        """

        vector = True if type(stations[0][3]) == list else False
        a_squared = 0.25 * (dx ** 2 + dy ** 2)
        lat = point[0]
        lon = point[1]
        elev = point[2]
        if topography is None:
            Ctag = C
        else:
            h = elev - topography
            Ctag = D / (-numpy.tanh((h - Hsl / 2) / b) / 2 + 0.5 + D / C)
        for n in range(len(stations)):
            if stations[n][0] == lat and stations[n][1] == lon and stations[n][2] == elev:
                return stations[n][3]
        wt = 0
        value = [0 for i in range(len(stations[0][3]))] if vector else 0
        for n in range(len(stations)):
            r = float((stations[n][0] - lat) ** 2 + (stations[n][1] - lon) ** 2) ** 0.5
            widw = 1.0 / (1. + (r ** 2) / a_squared)
            dz = math.fabs(stations[n][2] - elev)
            wedw = 1. / (1. + Ctag * (dz ** 2. / a_squared))
            wt += widw * wedw
        for n in range(len(stations)):
            r = float((stations[n][0] - lat) ** 2 + (stations[n][1] - lon) ** 2) ** 0.5
            widw = 1.0 / (1. + (r ** 2) / a_squared)
            dz = math.fabs(stations[n][2] - elev)
            wedw = 1. / (1. + Ctag * (dz ** 2. / a_squared))
            wb = widw * wedw
            if vector:
                for i in range(len(stations[0][3])):
                    value[i] += wb * stations[n][3][i] / wt
            else:
                value += wb * stations[n][3] / wt

        return value

    def checkInterpulation(self, stations, dx=20, dy=20, C=1000, D=5, Hsl=100, b=150):
        """
        Gets a list of stations, and performs an interpulation of the values in each station
        by the information of the other stations.
        Returns a list of the differences in percentage between the measured and interpulated values.
        """
        differences = []
        for station in stations:
            point = [station[0], station[1], station[2]]
            newStations = stations.copy()
            newStations.remove(station)
            inter = self.interp(point=point, topography=point[2], stations=newStations, dx=dx, dy=dy, C=C, D=D, Hsl=Hsl,
                                b=b)
            measure = station[3]
            difference = []
            for i in range(len(inter)):
                difference.append((inter[i] - measure[i]) / measure[i] * 100)
            differences.append(difference)

        return differences

    def interpPandas(self, points, stations, columnNames={"x": "x", "y": "y", "z": "z", "topography": "topography"},
                     dx=20, dy=20, C=1000, D=5, Hsl=100, b=150):

        points = points.reset_index()
        points["interpulation"] = None
        for i in range(len(points)):
            point = [points[columnNames["x"]][i], points[columnNames["y"]][i], points[columnNames["z"]][i]]
            topography = points[columnNames["topography"]][i] if columnNames["topography"] in points.columns else None
            points["interpulation"][i] = self.interp(point=point, stations=stations, topography=topography,
                                                     dx=dx, dy=dy, C=C, D=D, Hsl=Hsl, b=b)
        return points

    def interpArray(self, points, stations, columnNames={"x": "x", "y": "y", "z": "z", "topography": "topography"},
                    dx=20, dy=20, C=1000, D=5, Hsl=100, b=150):

        newPoints = pandas.DataFrame({"x": points[columnNames["x"]],
                                      "y": points[columnNames["y"]],
                                      "z": points[columnNames["z"]],
                                      "topography": points[columnNames["topography"]]})
        return self.interpPandas(points=newPoints, stations=stations, topography=newPoints["topography"],
                                 dx=dx, dy=dy, C=C, D=D, Hsl=Hsl, b=b)

    def windprofile(z, uref=3, href=24, he=24, lambdap=0.3, lambdaf=0.3, beta=0.3):

        # This function will return the wind velocity at height z

        # parameters:
        # z is the predicted height [m] above the ground
        # uref is the velocity reference [m/s] at height href
        # href is the reference height [m] of uref
        # he is the mean building height [m] above the ground
        # lambdap is the building percent when looking from top
        # lambdaf is the building percent when looking from the side
        # beta is the ratio between ustar and Uh

        # return
        # u is the predicted velocity [m/s] that is calculated in this function

        # variables:
        # l is the mixing length [m]

        ###############################
        ##    z=24
        #    uref=3
        #    href=24
        #    he=24
        #    lambdap=0.3
        #    lambdaf=0.3
        ##  assuming href is he
        uh = uref
        ###############################

        k = 0.41  # von Karman constant
        beta = beta  # For closed uniform natural canopies, Raupach 1996
        # beta = ustar / uh
        # LA commercial area has lambdap=0.28, lambdaf=0.27 he=24.5m, LC=66, Coceal 2004
        lc = he * (1 - lambdap) / lambdaf
        # calculating the mixing length at the canopy only to know d
        l = 2. * beta ** 3 * lc  # mixing length in the canopy
        d = l / k

        if href >= he:
            # l = k*(href-d) # mixing length above the canopy
            z01 = d / math.exp(k / beta)
            z0 = l / k * math.exp(-k / beta)  # from Omri
            #       print ('z0=',z0,z01,l, d/math.exp(k/beta))
            uh = uref * k / beta / math.log(((href - he) + d) / z0)
        else:  # href<he
            uh = uref / math.exp(beta * (href - he) / l)
            print(
                'there is no sense of giving uref below building mean height because we have to give mean velocity, above he, one point represent the mean value, below, it is not represent')

        #    print('he=',he)
        if z <= he:  # below the canopy
            #     l = 2. * beta**3 * lc # mixing length in the canopy
            u = uh * math.exp(beta * ((z - he)) / l)
        else:  # above the canopy
            #       calculate now for z
            #       l = k*(z-d) # mixing length above the canopy
            #       z0 = d/math.exp(k/beta)
            u = uh * beta / k * math.log(((z - he) + d) / z0)

        return u


if __name__ == "__main__":
    haifa1 = [
        [34.98535, 32.78562, 280],
        [34.99925, 32.77753, 370],
        [35.08025, 32.829758, 10],
        [34.99718, 32.81188, 75],
        [35.03992, 32.789, 18],
        [34.99066, 32.80078, 274],
        [35.00926, 32.78981, 150],
        [35.09444, 32.74263, 95],
        [34.97059, 32.79695, 215],
        [35.0554, 32.5119, 15],
        [35.02111, 32.791929, 104],
        [35.02026, 32.78685, 240],
        [35.04226, 32.76994, 90],
        [35.00167, 32.81644, 8],
        [35.03613, 32.73827, 507],
        [35.1123, 32.81148, 65],
        [35.077715, 32.814143, 25],
        [35.08511, 32.78866, 5],
        [35.1294, 32.7218, 201],
        [35.07873, 32.85182, 27],
        [35.091182, 32.854514, 15],
        [35.054671, 32.831263, 0],
        [35.00168, 32.78943, 190],
        [34.9651, 32.8226, 107]
    ]

    haifa0 = [8.7,
              6.6,
              35.7,
              23,
              16.8,
              10.6,
              10.9,
              8.8,
              4,
              8.2,
              13.6,
              8.4,
              10.6,
              50.7,
              5.8,
              8.3,
              20.6,
              8.7,
              10.9,
              19,
              19.9,
              24.6,
              8.7,
              5.7
              ]

    haifa0 = [7.4
        , 6.4
        , 12.3
        , 27.9
        , 12.7
        , 12.3
        , 6.6
        , 9.3
        , 6.6
        , 15.1
        , 9.8
        , 6.7
        , 35.4
        , 4.2
        , 4.5
        , 5.2
        , 4
        , 11.8
        , 7
        , 7.9
        , 6.1
        , 8.5
        , 10.5
              ]

    haifa0 = [6.6
        , 5.1
        , 13.2
        , 36.5
        , 9
        , 10.4
        , 5.8
        , 6
        , 5.9
        , 11.7
        , 6.6
        , 4.7
        , 46.5
        , 5.2
        , 5.8
        , 6.1
        , 3.2
        , 8.4
        , 7.8
        , 6.1
        , 6.8
        , 7.2
        , 11
              ]

    haifa0 = [19
        , 17.1
        , 33.7
        , 53.4
        , 14.8
        , 20.4
        , 21.1
        , 12.7
        , 19.7
        , 59
        , 32.9
        , 25
        , 17.6
        , 49
        , 8.5
        , 18.6
        , 23.6
        , 7.4
        , 17.3
        , 11.1
        , 21.5
        , 6.7
        , 26.5
        , 10.1
              ]

    # stations=[]
    # for i in range(len(haifa0)):
    #     stations.append([haifa1[i][0], haifa1[i][1], haifa0[i], haifa1[i][2]])
    # print (interp(35., 32.,  stations, elev = 111))
    # print (interp(stations[3][0], stations[3][1],  stations, elev = stations[3][2]))
