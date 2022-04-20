import math
import random

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

import numpy


class Particle:
    def __init__(self, x, v, m):
        self.x = x
        self.v = v
        self.m = m
        self.f = numpy.array([1, 1, 1])
        self.pot = 0
        self.eps = 120
        self.sigma = 3.38

    def kinenergy(self):
        return self.m * (self.v[0]**2 + self.v[1]**2 + self.v[2]**2) / 2

    def calculate_force(self, prt, n):
        self.f = numpy.array([0, 0, 0])
        self.pot = 0
        for i in range(len(prt)):
            if i != n:
                r = math.sqrt((self.x[0] - prt[i].x[0])**2 + (self.x[1] - prt[i].x[1])**2 + (self.x[2] - prt[i].x[2])**2)
                abf = 4 * 1.38 * self.eps * (-6 * self.sigma ** 6 / r ** 7 + 12 * self.sigma ** 12 / r ** 13)
                self.f[0] += abf * (self.x[0] - prt[i].x[0]) / r
                self.f[1] += abf * (self.x[1] - prt[i].x[1]) / r
                self.f[2] += abf * (self.x[2] - prt[i].x[2]) / r
                self.pot += 4 * 1.38 * self.eps * (-(self.sigma / r) ** 6 + (self.sigma / r) ** 12)


# class Spring:
#     def __init__(self, x, k):
#         self.x = x
#         self.k = k
#
#     def energy(self):
#         return self.k * (self.x[0]**2 + self.x[1]**2 + self.x[2]**2) / 2


# class Wall:
#     def __init__(self, a, x):
#         self.a = a
#         self.x = -x
#
#     def force(self, prt):
#         abx = math.sqrt(self.x[0]*self.x[0]+self.x[1]*self.x[1]+self.x[2]*self.x[2])
#         r = (self.x[0]*prt.x[0]+self.x[1]*prt.x[1]+self.x[2]*prt.x[2]+abx**2)/abx
#         abf = -6*self.a/r**7 + 12*self.a/r*13
#         f = [0, 0, 0]
#         f[0] = abf * self.x[0] / abx
#         f[1] = abf * self.x[1] / abx
#         f[2] = abf * self.x[2] / abx
#         return r, f


def integrate(prt, dt, tmax, xm):
    kin_t = [[]]
    poten = numpy.array([])
    pot = 0
    kinet = numpy.array([])
    kin = 0
    f = [[], []]
    tmin = 0
    steps = round((tmax - tmin) / dt)
    t = numpy.linspace(tmin, tmax, steps)
    coords = []
    for i in range(len(prt)):
        prt[i].calculate_force(prt, i)
        #f[i].append(math.sqrt(prt[i].f[0] ** 2 + prt[i].f[1] ** 2 + prt[i].f[2] ** 2))
        pot += prt[i].pot
        kin += prt[i].kinenergy()
    poten = numpy.append(poten, pot)
    kinet = numpy.append(kinet, kin)
    for i in range(len(prt)):
        coords.append([])
        coords[i].append(prt[i].x)
        coords[i].append(coords[i][0] + (tmax - tmin) / steps * prt[i].v)
        prt[i].x = coords[i][1]

    for i in range(len(prt)):
        prt[i].calculate_force(prt, i)
        #f[i].append(math.sqrt(prt[i].f[0] ** 2 + prt[i].f[1] ** 2 + prt[i].f[2] ** 2))
        pot += prt[i].pot
        kin += prt[i].kinenergy()
    poten = numpy.append(poten, pot)
    kinet = numpy.append(kinet, kin)
    t = [0]
    t.append((tmax - tmin) / steps)
    n_series = 20
    for i in range(steps - 2):
        if i > 1:
            kin = 0
            pot = 0
            for n in range(len(prt)):
                prt[n].calculate_force(prt, n)
                coords[n].append(2 * coords[n][i - 1] - coords[n][i - 2] + 0.8313*((tmax - tmin) / steps) ** 2 * (prt[n].f / prt[n].m))
                prt[n].x = coords[n][i]
                prt[n].v = (coords[n][i] - coords[n][i-1])/((tmax - tmin) / steps)
                kin += prt[n].kinenergy()
                if len(kin_t[len(kin_t)- 1]) < n_series:
                    kin_t[len(kin_t) - 1].append(prt[n].kinenergy())
                else:
                    kin_t.append([])
                pot += prt[n].pot
                #f[n].append(math.sqrt(prt[n].f[0] ** 2 + prt[n].f[1] ** 2 + prt[n].f[2] ** 2))
            poten = numpy.append(poten, pot)
            print(pot, kin)
            kinet = numpy.append(kinet, kin)
            t.append(t[i - 2] + (tmax - tmin) / steps)
    for i in range(steps - 2):
        for n in range(len(prt)):
            for j in range(3):
                if coords[n][i][j] != 0:
                    coords[n][i][j] = (coords[n][i][j] + xm) % (2 * xm) - xm
    return coords, t, f, poten, kinet, kin_t


def plot_update(data):
    ax.clear()
    xm = 15
    ax.set_xlim3d([-xm, xm])
    ax.set_ylim3d([-xm, xm])
    ax.set_zlim3d([-xm, xm])
    img = []
    for i in range(len(data)):
        img.append(ax.scatter3D(data[i][0], data[i][1], data[i][2]))
    return img


n = 20
prt = []
#for i in range(n):
#prt.append(Particle(numpy.array([0, 2, 1]), numpy.array([0, 0, 0]), 1))
part_mas = 39.9
xm = 15
T = 50
vm = 1.57*math.sqrt(T/part_mas)
for i in range(n):
    prt.append(Particle(numpy.array([random.uniform(-xm, xm), random.uniform(-xm, xm), random.uniform(-xm, xm)]), numpy.array([random.uniform(-vm, vm), random.uniform(-vm, vm), random.uniform(-vm, vm)]), part_mas))
data, t, f, poten, kinet, kin_t = integrate(prt, 0.001, 1, xm)
poten = poten / 2
coords = []
for i in range(len(data[1])):
    coords.append(math.sqrt((data[0][i][0] - data[1][i][0])**2 + (data[0][i][1] - data[1][i][1])**2 + (data[0][i][2] - data[1][i][2])**2))
plt.plot(t, kinet)
plt.plot(t, poten)
plt.plot(t, poten + kinet)
fig = plt.figure()
ax = p3.Axes3D(fig)
dataimg = []
step = []
for i in range(len(data[0])):
    if i % 10 == 0:
        for j in range(len(data)):
            step.append(data[j][i])
        dataimg.append(step)
    step = []
line_ani = animation.FuncAnimation(fig, plot_update, dataimg, interval = 100, blit=False)
plt.show()

print(len(kin_t))
while True:
    k = 0
    k = int(input())
    print(kin_t[k])
    n, bins, patches = plt.hist(kin_t[k], 10, density= True, facecolor='g')
    x = numpy.linspace(0, 1000, 1000)
    y = 4 * 3.14 * (1.66 * 10 ** -4 / (2 * 3.14 * 1.38 * T)) * x ** 2 * 94-math.exp(-1.66* 10 **-4 * x**2/(2*1.38*T))
    plt.plot()
    plt.show()


