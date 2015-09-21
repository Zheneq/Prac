import math
from rk import *


def f(t, y):
    return y


def f_(t, y):
    return y[0]


def f0(t, yv):
    return 2*yv[0] + math.exp(t) / math.cos(t) * (yv[1] + math.exp(t) * math.sin(t)**2)


def f1(t, yv):
    return yv[1] - 4. * math.cos(t) / math.exp(t) * yv[0]


def g0(t, yv):
    return yv[1] - 3*yv[0] + 9*t + 2


def g1(t, yv):
    return 6/7.*yv[0] - 2/7.*yv[1] + 4 * (t + 1)


h = 1e-3
y = 1.0
# # one equation
# for i in xrange(0, 100):
#     x = i * h
#     y = rungekutta(x, y, f, h)
#     print y - math.exp(x + h)

# # system of one equation
# y = [1.0]
# for i in xrange(0, 100):
#     x = i * h
#     y = rungekuttavec(x, y, [f_], h)
#     print y[0] - math.exp(x + h)

# system of two equations
fv = [f0, f1]
yv = [0., 1.]
for i in xrange(1, 100):
    x = i * h
    y = rungekuttavec(x, yv, fv, h)
    print y[0] - math.exp(2*x) * math.sin(x)
    print y[1] - math.exp(x) * math.cos(2*x)

# # system of two equations
# fv = [g0, g1]
# yv = [-1., 4.]
# for i in xrange(1, 100):
#     x = i * h
#     y = rungekuttavec(x, yv, fv, h)
#     print y[0] - x**2 - 3*x + 1
#     print y[1] - 3*x**2 - 2*x - 4

