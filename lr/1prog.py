# Задание:
# r(t) = 2 + cos(6t)
# phi(t) = t + sin(6t)
# найти траекторию, поставить векторочки центростремительного ускорения и скорости

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as l2
from matplotlib.animation import FuncAnimation
import sympy as sp
import math

def rotate_2d (X, Y, Alpha):  #Уравнения, определяющие преобразование в двух измерениях, которое поворачивает оси xy против часовой стрелки под углом theta  в оси x'y'
    RX = X*np.cos(Alpha) - Y*np.sin(Alpha)
    RY = Y*np.cos(Alpha) + X*np.sin(Alpha)
    return RX, RY

#ПЕРЕМЕННЫЕ И КОНСТАНТЫ
#понятные
frame_frequency = 10 # частота кадров
frame_num = 1000 #количество кадров?
t = sp.Symbol('t') # символ тэ

#непонятные
v_scl = 0.1
a_scl = v_scl**2
ar_scl = 0.025
t_end = 10 #*math.pi/6 #10*pi/6 - конец "времени"

#данные значения
r = 2 + sp.cos(6*t)
phi = t + sp.sin(6*t)

#переход к полярным координатам
x = r * sp.sin(phi)
y = r * sp.cos(phi)
#скорость как произвдная по координате
v_x = sp.diff(x, t)
v_y = sp.diff(y, t)
#ускорение как производная скорости
a_x = sp.diff(v_x, t)
a_y = sp.diff(v_y, t)

#списки значений
T = np.linspace(0, t_end, frame_num) #(от, до, num - сколько чисел среди этих значений) //интервал дедить на num равных частей numpy.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None)
X = np.zeros_like(T) #Функция zeros_like() возвращает новый массив из нулей с формой и типом данных указанного массива T.T - существующий массив
Y = np.zeros_like(T)
V_X = np.zeros_like(T)
V_Y = np.zeros_like(T)
A_X = np.zeros_like(T)
A_Y = np.zeros_like(T)
R_X = np.zeros_like(T)
R_Y = np.zeros_like(T)
R_ = np.zeros_like(T)

# заполнения их значениями по формулам
for i in np.arange(len(T)):
    X[i] = sp.Subs(x, t, T[i]) #(берем выражение x, и везде где видим в нем t, заменяем его на T[i]) и так дня каждого значения по массиву из 1000 элементов(кадров) (можно еще больше значений (expression(x,y), (x,y), (3,4))- x,y заменятся на 3,4)
    Y[i] = sp.Subs(y, t, T[i])
    V_X[i] = sp.Subs(v_x, t, T[i])
    V_Y[i] = sp.Subs(v_y, t, T[i])
    A_X[i] = sp.Subs(a_x, t, T[i])
    A_Y[i] = sp.Subs(a_y, t, T[i])

# ничего непонятно в следующием абзаце
    R_X[i] = (V_X[i] ** 2 + V_Y[i] ** 2) / (A_X[i] * V_Y[i] - A_Y[i] * V_X[i])  # хз что это
    R_Y[i] = R_X[i] * V_Y[i]  # это вот просто переделывание на Y из того, что выбо выше
    R_X[i] *= V_X[i]  # блин, какая-то лажа
    R_[i] = math.sqrt(R_X[i] ** 2 + R_Y[i] ** 2)

V_X = V_X * v_scl
V_Y = V_Y * v_scl
A_X = A_X * a_scl
A_Y = A_Y * a_scl

fig = plt.figure() # создание одласти для рисования
fig.set_facecolor('#c1c6fc') #цвет снаружи поля

ax1 = fig.add_subplot(1, 1, 1)  # задаём кол-во участков для рисования
ax1.set_facecolor('#FFF0F5') #цвет внутри осей
ax1.axis('equal')
ax1.set(xlim=[-3, 3], ylim=[-3, 3])

# ax1.plot(X, Y)
# ax1.plot([X.min(), X.max()], [0, 0], 'black') #хз что это. явно что-то черного цвета

P, = ax1.plot(X[0], Y[0], marker = 'o', color = 'xkcd:green blue') #o - точка; синтаксис: plot(x, y, цвет и вид того, что отображаем)
V, = ax1.plot([X[0], X[0] + V_X[0]], [Y[0], Y[0] + V_Y[0]], 'xkcd:cornflower') # палка вектора скорости
R, = ax1.plot(X[0], Y[0], color = 'xkcd:dark mint') #кажется будто бы это лишняя палка посередине. ну да.  м. нет. без нее вообще ничего не работает
A, = ax1.plot([X[0], X[0] + A_X[0]], [Y[0], Y[0] + A_Y[0]], 'xkcd:carnation pink') # зеленая палка вектора ускорения

######################################################
#Цвета xkcd: https://xkcd.com/color/rgb/ и https://russianblogs.com/article/69592450442/
#linestyles: https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html

# Alpha = np.linspace(0, math.pi*2, 100)
# Circ, = ax1.plot(X[0] + R_Y[0] - R_[0]*np.cos(Alpha),
#                 Y[0] - R_X[0] - R_[0]*np.sin(Alpha),
#                 color = 'xkcd:mint', linestyle = '-.')
#здесь почему-то остается начальное положение круга. интересно.

#А здесь уже не остается
Alpha = np.linspace(0, 6.28, 100)
Circ, = ax1.plot(X[0]+R_[0]*R_X[0] * np.cos(Alpha), Y[0]+R_[0]*R_Y[0] * np.sin(Alpha), color = '#7B68EE', linestyle = '-.')

ArrowX = np.array([-2*ar_scl, 0, -2*ar_scl]) #что-то важное без чего ничего не работеат
ArrowY = np.array([ar_scl, 0, -ar_scl]) #и это тоже

VArrowX, VArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(V_Y[0], V_X[0]))
VArrow, = ax1.plot(VArrowX + X[0] + V_X[0], VArrowY + Y[0] + V_Y[0], 'xkcd:cornflower') #наконечник вектора скорости

RArrowX, RArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(Y[0], X[0]))
RArrow, = ax1.plot(RArrowX + X[0], RArrowY + Y[0], color = 'black')

AArrowX, AArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(A_Y[0], A_X[0]))
AArrow, = ax1.plot(AArrowX + X[0] + A_X[0], AArrowY + Y[0] + A_Y[0], 'xkcd:carnation pink') # наконечник вектора ускорения

ax1.plot(X, Y, 'xkcd:soft blue', linestyle = 'solid') #траектория
def anim (i):
    P.set_data(X[i], Y[i])
    V.set_data([X[i], X[i] + V_X[i]],
               [Y[i], Y[i] + V_Y[i]])
    R.set_data([0, X[i]], [0, Y[i]])
    A.set_data([X[i], X[i] + A_X[i]],
               [Y[i], Y[i] + A_Y[i]])
    VArrowX, VArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(V_Y[i], V_X[i]))
    VArrow.set_data(VArrowX + X[i] + V_X[i],
                    VArrowY + Y[i] + V_Y[i])
    RArrowX, RArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(Y[i], X[i]))
    RArrow.set_data(RArrowX + X[i], RArrowY + Y[i])
    AArrowX, AArrowY = rotate_2d(ArrowX, ArrowY, math.atan2(A_Y[i], A_X[i]))
    AArrow.set_data(AArrowX + X[i] + A_X[i],
                    AArrowY + Y[i] + A_Y[i])
    Circ.set_data(X[i] + R_Y[i] - R_[i]*np.cos(Alpha),
                  Y[i] - R_X[i] - R_[i]*np.sin(Alpha))
    return P, V, VArrow, R, RArrow, A, AArrow, Circ,

anim1 = FuncAnimation(fig, anim, frames = frame_num, interval = frame_num/frame_frequency, blit = True)

custom_lines = [l2([0], [0], color='xkcd:soft blue', linestyle='solid'),
                l2([0], [0], color='xkcd:dark mint', lw=4),
                l2([0], [0], color='xkcd:cornflower', lw=4),
                l2([0], [0], color='xkcd:carnation pink', lw=4),
                l2([0], [0], color='#7B68EE', linestyle='-.', lw=2)]
ax1.legend(custom_lines, ['Траектория', 'Радиус-вектор', 'Вектор скорости', 'Ускорение', 'Кружок:3'], loc='lower left')

plt.show()
