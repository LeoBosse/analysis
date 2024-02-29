#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import observation as obs
from geometry import GetLonLatFromName
from itertools import combinations

DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371 # #km

###########
# For a list of polarisation observations (AoLPs and instrument setup), computes the theoretical planes of the current (perpendicular or parallel to the AoLP) for each one, and find their intersections. 
# Prints the intersections, and average intersection with the angle difference between each intersection.
# You want the difference to be as close to zero as possible.
###########


# Fill in a list of dictionnaries with the AoLP and instrumental setup for each observations used in this triangulation
gdcu = {
'name': "lacsud",
'h': 110, #km
'a': 12 * DtoR,
'e': 75 * DtoR,
'AoLP': 90 * DtoR #90
}

corbel = {
'name': "skibotn",
'h': 110, #km
'a': 282 * DtoR,
'e': 44 * DtoR,
'AoLP': 0  * DtoR #-8
}

# carmen = {
# 'name': "skibotn",
# 'h': 110, #km
# 'a': -97 * DtoR,
# 'e': 45 * DtoR,
# 'AoLP': 40 * DtoR
# }

#List of the dictionnaries to use
instruments = [gdcu, corbel]


#Complete the dictionnaries with the observation objects and perp and para planes
for i in instruments:
    i['lon'] = GetLonLatFromName(i['name'])[0]
    i['lat'] = GetLonLatFromName(i['name'])[1]
    i['observation'] = obs.ObservationPoint(i['lon'], i['lat'], i['h'], i['a'], i['e'])
    i['Pperp'] = i['observation'].GetPolaPlane(i['AoLP'], 'perp', 'A') #Norm of the plane that is perp to the los and the AoLP (cross product)
    i['Ppara'] = i['observation'].GetPolaPlane(i['AoLP'], 'para', 'A') #Norm of the plane that is parallel to the AoLP (contains AoLP and los -> norm is the cross product)


def UENToAzEl(uen):
    '''Converts up-east-north coordinates to azimuth and elevation angles.'''
    u, e, n = uen
    az = np.arctan2(e, n)
    el = np.arcsin(u)
    return az, el


# All possible pairings of the instruments. (1, 2)==(2, 1)
pairings = combinations(instruments, 2)

perp_intersections_uen = []
perp_intersections_azel = []
para_intersections_uen = []
para_intersections_azel = []

# Compute the intersection of the plane for each pair of observation
for i, pair in enumerate(pairings):
    perp_uen = np.cross(pair[0]['Pperp'], pair[1]['Pperp'])
    perp_uen /= np.linalg.norm(perp_uen)

    para_uen = np.cross(pair[0]['Ppara'], pair[1]['Ppara'])
    para_uen /= np.linalg.norm(para_uen)

    perp_azel = UENToAzEl(perp_uen)
    para_azel = UENToAzEl(para_uen)

    perp_intersections_uen.append(perp_uen)
    para_intersections_uen.append(para_uen)
    perp_intersections_azel.append(perp_azel)
    para_intersections_azel.append(para_azel)


if len(instruments) > 2:
    perp_diff = []
    para_diff = []
    perp_pairings = combinations(perp_intersections_uen, 2)
    para_pairings = combinations(para_intersections_uen, 2)

    print('Angles between every intersections of the perpendicular planes')
    for i, pair in enumerate(perp_pairings):
        angle_diff = np.arccos(np.dot(pair[0], pair[1]))
        perp_diff.append(angle_diff)
        print(angle_diff*RtoD)

    print('Angles between every intersections of the parallel planes')
    for i, pair in enumerate(para_pairings):
        angle_diff = np.arccos(np.dot(pair[0], pair[1]))
        para_diff.append(angle_diff)
        print(angle_diff*RtoD)

    perp_error = np.average(perp_diff)
    para_error = np.average(para_diff)
    if perp_error > para_error:
        print(f'PARALLEL currents are closer, with an average difference angle of {para_error*RtoD}')
    elif perp_error < para_error:
        print(f'PERPENDICULAR currents are closer, with an average difference angle of {perp_error*RtoD}')
    else:
        print('Parallel and perpendicular have exactly the same average difference! Hard to tell which is better...')

    print('All intersections for perp current')
    # print(perp_intersections_uen)
    for i in range(len(perp_intersections_uen)):
        print(perp_intersections_uen[i])

    print('All intersections for para current')
    # print(para_intersections_uen)
    for i in range(len(para_intersections_uen)):
        print(para_intersections_uen[i])


perp_current = np.average(perp_intersections_uen, axis=0)
perp_current /= np.linalg.norm(perp_current)
print("Average perpendicular current (uen)", perp_current)
print("Average perpendicular current (azel)", np.array(UENToAzEl(perp_current))*RtoD)

para_current = np.average(para_intersections_uen, axis=0)
para_current /= np.linalg.norm(para_current)
print("Average parallel current (uen)", para_current)
print("Average parallel current (azel)", np.array(UENToAzEl(para_current))*RtoD)
