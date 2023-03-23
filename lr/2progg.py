# Задание:
# r(t) = 2 + cos(6t)
# phi(t) = t + sin(6t)
# найти траекторию, поставить векторочки центростремительного ускорения и скорости

import copy # импортируем модуль copy для создания копий объектов
import math # импортируем модуль math для математических функций
import random # импортируем модуль random для генерации случайных чисел

class Point: # определяем класс Point для представления точек на плоскости
    def __init__(self, x=0, y=0): # конструктор класса, принимает координаты x и y
        self.x = x # сохраняем координату x в атрибуте self.x
        self.y = y # сохраняем координату y в атрибуте self.y

    def __str__(self): # метод для преобразования объекта в строку
        return f"({self.x}, {self.y})" # возвращаем строку в формате (x, y)

    def __add__(self, other): # метод для сложения двух точек
        return Point(self.x + other.x, self.y + other.y) # возвращаем новую точку с суммой координат

    def __sub__(self, other): # метод для вычитания двух точек
        return Point(self.x - other.x, self.y - other.y) # возвращаем новую точку с разностью координат

    def __mul__(self, other): # метод для умножения точки на число или на другую точку
        if isinstance(other, (int, float)): # если other - число
            return Point(self.x * other, self.y * other) # возвращаем новую точку с произведением координат на число
        elif isinstance(other, Point): # если other - точка
            return self.x * other.x + self.y * other.y # возвращаем скалярное произведение двух точек
        else: # иначе
            raise TypeError("unsupported operand type(s) for *") # выбрасываем исключение о неподдерживаемом типе операнда

    def __rmul__(self, other): # метод для умножения числа на точку (для коммутативности)
        return self * other # возвращаем результат умножения точки на число

    def __truediv__(self, other): # метод для деления точки на число
        if isinstance(other, (int, float)): # если other - число
            return Point(self.x / other, self.y / other) # возвращаем новую точку с частным координат на число
        else: # иначе
            raise TypeError("unsupported operand type(s) for /") # выбрасываем исключение о неподдерживаемом типе операнда

    def length(self): # метод для вычисления длины вектора от начала координат до точки
        return math.sqrt(self * self) # возвращаем корень из скалярного произведения точки на себя

    def angle(self): # метод для вычисления угла между вектором от начала координат до точки и осью OX
        return math.atan2(self.y, self.x) # возвращаем арктангенс отношения координаты y к x

    def rotate(self, angle): # метод для поворота точки на заданный угол относительно начала координат
        cos = math.cos(angle) # вычисляем косинус угла
        sin = math.sin(angle) # вычисляем синус угла
        return Point(self.x * cos - self.y * sin, self.x * sin + self.y * cos) # возвращаем новую точку с координатами после поворота

    def polar(self): # метод для перевода точки из декартовых координат в полярные
        return (self.length(), self.angle()) # возвращаем пару (расстояние, угол)

    def cartesian(self): # метод для перевода точки из полярных координат в декартовые
        return Point(self.x * math.cos(self.y), self.x * math.sin(self.y)) # возвращаем новую точку с координатами x = r * cos(phi), y = r * sin(phi)

class Curve: # определяем класс Curve для представления кривых на плоскости
    def __init__(self, r, phi): # конструктор класса, принимает функции r(t) и phi(t)
        self.r = r # сохраняем функцию r(t) в атрибуте self.r
        self.phi = phi # сохраняем функцию phi(t) в атрибуте self.phi

    def __call__(self, t): # метод для вычисления точки на кривой по параметру t
        return Point(self.r(t), self.phi(t)).cartesian() # возвращаем точку с полярными координатами (r(t), phi(t)), переведенную в декартовы

    def points(self, n): # метод для генерации n равномерно распределенных точек на кривой
        for i in range(n): # для каждого i от 0 до n-1
            t = i / (n - 1) * 2 * math.pi # вычисляем параметр t от 0 до 2*pi
            yield self(t) # возвращаем точку на кривой с параметром t

    def plot(self): # метод для рисования кривой на графике
        import matplotlib.pyplot as plt # импортируем модуль matplotlib.pyplot для работы с графикой
        x = [] # создаем пустой список для координат x
        y = [] # создаем пустой список для координат y
        for p in self.points(100): # для каждой точки из 100 точек на кривой
            x.append(p.x) # добавляем координату x в список x
            y.append(p.y) # добавляем координату y в список y
        plt.plot(x, y) # рисуем график по спискам x и y
        plt.axis("equal") # делаем оси равномерными по масштабу
        plt.show() # показываем график

    def velocity(self, t): # метод для вычисления скорости точки на кривой по параметру t
        h = 0.0001 # задаем маленькое число h для приближенного дифференцирования
        return (self(t + h) - self(t)) / h # возвращаем приближенное значение производной по t

    def acceleration(self, t): # метод для вычисления ускорения точки на кривой по параметру t
        h = 0.0001 # задаем маленькое число h для приближенного дифференцирования
        return (self.velocity(t + h) - self.velocity(t)) / h # возвращаем приближенное значение производной по t скорости

curve = Curve(lambda t: 2 + math.cos(6*t), lambda t: t + math.sin(6*t)) 
# создаем объект класса Curve с функциями r(t) = 2 + cos(6t), phi(t) = t + sin(6t)

curve.plot() 
# рисуем кривую на граф
