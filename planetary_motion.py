import numpy as np
import datetime
import math
from scipy import optimize
# (perihelion, period)
orbital_data = {'earth':({2018:(1,3), 2019:(1,3),2020:(1,5)},
                365.25635, 0.016710219,1.496e8)}

# (Axial tilt, rotation (rad/day) per day)
rotation_data = {'earth':(23.44, 6.30046407)}

'''Orbital up vector (orthogonal to earths orbit plane)
Note I am using plane span({(1,0,0), (0,0,1)})
as the plane of reference for all computations. 
'''
orbital_up = np.array([0,1,0])

def mean_anomaly(date:"(Year, Month, Day)", body = 'earth'):
    iyear, imonth, iday = date
    perihelion_data, period, eccentricity, semi_major_axis = orbital_data[body]
    perihelion = perihelion_data[iyear]

    day = datetime.date(iyear, imonth, iday)
    perihelion_day = datetime.date(iyear, perihelion[0], perihelion[1])

    #If the date is before the this years perihelion date, 
    #use last years perihelion
    if(day < perihelion_day):
        perihelion = perihelion_data[iyear - 1]
        perihelion_day = datetime.date(iyear - 1, perihelion[0], perihelion[1])

    t = (day - perihelion_day).days
    n = 2*(math.pi) / period
    return n * t

def eccentric_anomaly(date:"(Year, Month, Day)", body = 'earth'):
    perihelion_data, period, eccentricity, semi_major_axis = orbital_data[body]
    M = mean_anomaly(date, body)

    keplerf = lambda E: E - eccentricity*math.sin(E) - M
    keplerf_prime = lambda E: 1 -eccentricity*math.cos(E)
    keplerf_prime2 = lambda E: eccentricity*math.sin(E)

    E = optimize.newton(keplerf, 0,
                        fprime = keplerf_prime,
                        fprime2 = keplerf_prime2)
    return E

def true_anomaly(date:"(Year, Month, Day)", body = 'earth'):
    perihelion_data, period, eccentricity, semi_major_axis = orbital_data[body]
    E = eccentric_anomaly(date, body)

    Theta = 2 * math.atan(math.sqrt((1+eccentricity)/(1-eccentricity)) *
            math.tan(E/2))
    return Theta

def heliocentric_distance(date:"(Year, Month, Day)", body = 'earth'):
    perihelion_data, period, eccentricity, semi_major_axis = orbital_data[body]

    E = eccentric_anomaly(date, body)
    Theta = true_anomaly(date,body)
    r = semi_major_axis * (1 - eccentricity*math.cos(E))
    return r

def cartesian_coordinates(date:"(Year, Month, Day)", body = 'earth'):
    Theta = true_anomaly(date,body)
    r = heliocentric_distance(date,body)
    return np.array([r * math.cos(Theta),0,r * math.sin(Theta)])
