import numpy as np
from math import cos, sin
from scipy.integrate import odeint

def solve_system(y, t, m1, m2, r, J, l, alpha, a, g):
	phi = y[0]
	thetta = y[1]
	dy = np.zeros(4)

	dy[0] = y[2]
	dy[1] = y[3]

	a11 = m1 * (a ** 2 + r ** 2 - 2 * a * r * cos(phi) + J  + m2 * r ** 2)
	a12 = m2 * l * r * cos(alpha + thetta)
	a21 = r * cos(alpha + thetta)
	a22 = l

	b1 = -m1 * a * r * sin(phi) * y[2] ** 2
	b1 += m2 * l * r * sin(alpha + thetta) * y[3] ** 2
	b1 += (m1 + m2) * g * r * sin(alpha) - m1 * g * a * sin(alpha + phi)
	
	b2 = -g * sin(thetta)

	det_a = a11 * a22 - a21 * a12
	dy[2] = (b1 * a22 - b2 * a12) / det_a
	dy[3] = (a11 * b2 - a21 * b1) / det_a

	return dy

def get_coords(y0, t, m1, m2, r, J, l, alpha, a, g):
	Y = odeint(solve_system, y0, t, (m1, m2, r, J, l, alpha, a, g))

	_phi = Y[:, 0]
	_psi = Y[:, 1]

	return _phi, _psi

def time_mapper(times, values):
	_map = dict(zip(times, values))

	def mapper(t):
		return _map[t]

	return mapper

