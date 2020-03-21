from scipy import misc
import numpy as np
from numpy.linalg import inv

def rotation_ecef2enu(theta_A,phi_A):
	# entrie: 
	# - position theta_A, phi_A (in radian) of the point A on the sphere
	# outputs: 
	# - rotation matrix M_ecef2A from ECEF frame to (East, North, Up) frame "enu" in point A
	# NB: inverse operation (enu 2 ecef) = transpose !!
	# Laundal & Richmond, SSR, 2016

	M_ecef2enuA = np.zeros((3,3)) 
	M_ecef2enuA[0,0] = -np.sin(phi_A)
	M_ecef2enuA[0,1] = np.cos(phi_A)
	M_ecef2enuA[0,2] = 0.0
	M_ecef2enuA[1,0] = -np.cos(theta_A)*np.cos(phi_A)
	M_ecef2enuA[1,1] = -np.cos(theta_A)*np.sin(phi_A)
	M_ecef2enuA[1,2] = np.sin(theta_A)
	M_ecef2enuA[2,0] = np.sin(theta_A)*np.cos(phi_A)
	M_ecef2enuA[2,1] = np.sin(theta_A)*np.sin(phi_A)
	M_ecef2enuA[2,2] = np.cos(theta_A)

	return M_ecef2enuA

def dot_product(u,v):
	# entrie: 
	# - vectors u and v
	# outputs: 
	# - dot product dot_prod = u.v

	dot_prod = u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

	return dot_prod

def cross_product(u,v):
	# entrie: 
	# - vectors u and v
	# outputs: 
	# - cross product cross_prod = uxv

	cross_prod = np.zeros((3,1))
	cross_prod[0] = u[1]*v[2]-u[2]*v[1]
	cross_prod[1] = u[2]*v[0]-u[0]*v[2]
	cross_prod[2] = u[0]*v[1]-u[1]*v[0]

	return cross_prod

