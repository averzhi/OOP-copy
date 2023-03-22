import numpy as np
from math import pi, cos, sin
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
from dynamics import time_mapper, solve_system

# Параметры для кадров

BEGIN_TIME = 0
END_TIME = 10
NUMBER_OF_FRAMES = 1001
FULL_CIRCLE = 360

# -------ФУНКЦИИ-------

# Поворачивает вектор относительно нулевой точки координат
def rotate_vector(vector, angle):
	return np.matmul(np.array([[cos(angle), -sin(angle)],[sin(angle), cos(angle)]]), vector)

# Принимает на вход матрицу, у которой каждая строка - вектор, и разворачивает все вектры в этом массиве
def rotate_vectors_in_array(vector_matrix, angle):
	return np.array([rotate_vector(vector_matrix[:, i], angle) for i in range(vector_matrix.shape[1])]).transpose()

# Просчитывание точек с помощью ф-ции rule. Аргументы из values
def dots_amount(rule, values):
	return [rule(value) for value in values]

# -------КОНСТАНТЫ(отображение)-------
plane_length = 3 # длина плоскости
plane_first_position = np.array([[4], [1]]) # точка начала координат

disk_first_position = -3 # начальное положение диска
number_of_spokes = 10 # количество спиц

# -------КОНСТАНТЫ-------
m1 =1 # масса диска
m2 = 0.1 # масса груза
r = 0.5 # радиус диска
J = 0.2 # момент инерции диска относительно центра масс
l = 1 # длина стержня
alpha = pi / 12  # угол между наклонной плоскостью и горзионтом
a = 0.25 # расстояние от центра масс до геометрического центра диска
g = 9.8 # ускорение свободного падения

y0 = [0, 0, -pi, pi] # начальные условия

#-----------------------------------------------

# -------ЗАКОНЫ-------
t = np.linspace(0, END_TIME, NUMBER_OF_FRAMES)
model_consts = (m1, m2, r, J, l, alpha, a, g)
y = odeint(solve_system, y0, t, model_consts)
phi, psi, dphi, dpsi = map(lambda index: time_mapper(t, y[:, index]), range(4))
ddphi = time_mapper(t, [solve_system(y0, time, *model_consts)[2] for y0, time in zip(y, t)])
ddpsi = time_mapper(t, [solve_system(y0, time, *model_consts)[3] for y0, time in zip(y, t)])

Rx = lambda t: (m1 + m2) * r * ddphi(t) - m1 * a * (ddphi(t) * cos(phi(t)) - dphi(t) ** 2 * sin(phi(t))) + m2 * l * (ddpsi(t) * cos(alpha + psi(t)) - dpsi(t) ** 2 * sin(alpha + psi(t))) + (m1 + m2) * g * sin(alpha)
Ry = lambda t: m1 * a * (ddphi(t) + sin(phi(t)) + dphi(t) ** 2 * cos(phi(t))) + m2 * l * (ddphi(t) * sin(alpha + psi(t)) + dpsi(t) ** 2 * cos(alpha + psi(t))) + (m1 + m2) * g * cos(alpha)

# -------ФУНКЦИЯ ДЛЯ ГЕНЕРИРОВАНИЯ ТОЧЕК ДЛЯ ОБЬЕКТОВ-------

# Конечные точки плоскости
def rule_for_inclined_plane_dots(t):
	return plane_first_position + np.array([[-plane_length, 0, -plane_length * cos(alpha)], [0, 0, plane_length * sin(alpha)]])

# Центр диска
def rule_for_disk_center(t):
	disk_center = np.array([[disk_first_position + r * phi(t)], [r]]) # локальные координаты
	return disk_center

# Центр массы диска
def rule_for_mass_center(t):
	mass_center = np.array([[0], [-a]])
	mass_center = rotate_vector(mass_center, -phi(t)) + rule_for_disk_center(t) # локальные координаты
	mass_center = rotate_vector(mass_center, -alpha) + plane_first_position # глобальные координат
	return mass_center

# Точки спиц окружности
def rule_for_spokes_dots(t):
	spokes_angles = np.linspace(0, 2 * pi, number_of_spokes)
	spokes_dots = np.zeros([2, 2 * number_of_spokes])
	for i, a in enumerate(spokes_angles):
		spokes_dots[0, 2 * i] = r * np.cos(a)
		spokes_dots[1, 2 * i] = r * np.sin(a)
	spokes_dots = rotate_vectors_in_array(spokes_dots, -phi(t))
	spokes_dots += rule_for_disk_center(t) # локальные координаты
	spokes_dots = rotate_vectors_in_array(spokes_dots, -alpha) + plane_first_position # глобальные координаты
	return spokes_dots

# Точки стержня
def rule_for_kernel_dots(t):
	kernel_dots = np.array([[0, 0], [0, -l]])
	kernel_dots = rotate_vectors_in_array(kernel_dots, alpha + psi(t)) # локальные координаты относительно центра
	kernel_dots += rule_for_disk_center(t) # локальные координаты
	kernel_dots = rotate_vectors_in_array(kernel_dots, -alpha) + plane_first_position # глобальные координаты
	return kernel_dots

# Точки окружности диска
def rule_for_disk_dots(t):
	disk_angles = np.linspace(0, 2 * pi, FULL_CIRCLE)
	disk_dots = r * np.array([np.cos(disk_angles), np.sin(disk_angles)]) # локальные координаты относительно центра
	disk_dots += rule_for_disk_center(t) # локальные координаты
	disk_dots = rotate_vectors_in_array(disk_dots, -alpha) + plane_first_position # глобальные координаты
	return disk_dots

# ------ГОТОВЫЕ ТОЧКИ ГРАФ МОДЕЛЕЙ-------

inclined_plane_dots = dots_amount(rule_for_inclined_plane_dots, t)
disk_dots = dots_amount(rule_for_disk_dots, t)
mass_center_dots = dots_amount(rule_for_mass_center, t)
kernel_dots = dots_amount(rule_for_kernel_dots, t)
spokes_dots = dots_amount(rule_for_spokes_dots, t)

#-----------------------------------------------

# Отрисовка графиков (рендеринг)

graph_figure = plt.figure(figsize=[2, 4])

# Добавляем оси
phi_graph = graph_figure.add_subplot(2, 4, 1, title="$\\varphi(t)$") # title - значки над осями
dphi_graph = graph_figure.add_subplot(2, 4, 2, title="$\\dot{\\varphi}(t)$") # 2 строки, 4 столбца, 2е место
ddphi_graph = graph_figure.add_subplot(2, 4, 3, title="$\\ddot{\\varphi}(t)$")

psi_graph = graph_figure.add_subplot(2, 4, 4, title="$\\psi(t)$")
dpsi_graph = graph_figure.add_subplot(2, 4, 5, title="$\\dot{\\psi}(t)$")
ddpsi_graph = graph_figure.add_subplot(2, 4, 6, title="$\\ddot{\\psi}(t)$")

rx_graph = graph_figure.add_subplot(2, 4, 7, title="$R_x$")
ry_graph = graph_figure.add_subplot(2, 4, 8, title="$R_y$")

# Рисуем сами графики функций

phi_graph.plot(t, dots_amount(phi, t))
dphi_graph.plot(t, dots_amount(dphi, t))
ddphi_graph.plot(t, dots_amount(ddphi, t))

psi_graph.plot(t, dots_amount(psi, t))
dpsi_graph.plot(t, dots_amount(dpsi, t))
ddpsi_graph.plot(t, dots_amount(ddpsi, t))

rx_graph.plot(t, dots_amount(Rx, t))
ry_graph.plot(t, dots_amount(Ry, t))

# -------ГРАФИЧЕСКИЕ ПРЕДСТАВЛЕНИЯ-------

figure = plt.figure(figsize=[4, 4])
subplot = figure.add_subplot(1, 1, 1)
subplot.axis("equal")
subplot.set(xlim=[0, 5], ylim=[0, 5])

inclined_plane, = subplot.plot(*inclined_plane_dots[0], color="#faee66")
disk, = subplot.plot(*disk_dots[0], color="#8e82fe")
spokes, = subplot.plot(*spokes_dots[0], color="#8e82fe")
mass_center, = subplot.plot(*mass_center_dots[0], color="#ff5b00", marker="o")
kernel, = subplot.plot(*kernel_dots[0], color="#fe01b1")

#-----------------------------------------------

# -------ФУНКЦИЯ АНИМАЦИИ-------

def animate_rule(i):
	kernel.set_data(kernel_dots[i])
	disk.set_data(disk_dots[i])
	inclined_plane.set_data(inclined_plane_dots[i])
	spokes.set_data(spokes_dots[i])
	mass_center.set_data(mass_center_dots[i])

	return kernel, disk, inclined_plane

animation = FuncAnimation(figure, animate_rule, frames=NUMBER_OF_FRAMES, interval=60)
plt.show()
