#!-*-coding:utf-8-*-
__author__ = 'Root'

from tdma import *
import math

PROB_DEF = 0
PROB_LOC = 1

SCH_LIN = 0
SCH_NONLIN = 1
SCH_ITER = 2

prob_def_y0 = lambda x: math.exp(-(x - 5.0)**2)
prob_def_k  = lambda t, x, y: 1.0
prob_def_f  = lambda t, x, y: 0.0
prob_def_m1 = lambda t: 0.0
prob_def_m2 = lambda t: 0.0

prob_loc_y0 = lambda x: 0.0
prob_loc_k  = lambda t, x, y, k0=1.0, s=3.0: k0 * y**s
prob_loc_f  = lambda t, x, y: 0.0
prob_loc_m1 = lambda t, a0=2.0, tf=100.0, n=-1./3.: a0 * (tf - t)**n
prob_loc_m2 = lambda t: 0.0

probs = [(prob_def_y0, prob_def_k, prob_def_f, prob_def_m1, prob_def_m2),
         (prob_loc_y0, prob_loc_k, prob_loc_f, prob_loc_m1, prob_loc_m2)]


def buildMatrix_linear(nx, h, tau, t, sigma, y, k, f):
    _a = [sigma / h**2] * nx
    _b = [sigma / h**2] * nx
    _c = [-((1 / tau) + (2 * sigma / h**2))] * nx
    _f = [-f(t, ix * h, y[ix]) - (1 - sigma) * (y[ix+1] - 2 * y[ix] + y[ix-1]) / h**2 - y[ix] / tau for ix in range(nx)]
    return _a, _b, _c, _f


def buildMatrix_nonlinear(nx, h, tau, t, sigma, y, k, f):
    # _a = [a(t, ix * h, y[ix]) * sigma / h**2 for ix in range(nx)]
    # _b = [a(t, (ix+1) * h, y[ix+1]) * sigma / h**2 for ix in range(nx)]
    # _c = [-((rho(t, ix * h, y[ix]) / tau) + ((a(t, ix * h, y[ix]) + a(t, (ix+1) * h, y[ix+1])) * sigma / h**2)) for ix in range(nx)]
    # _f = [-f(t, ix * h, y[ix]) -
    #       ((1 - sigma) / h**2) * (a(t, (ix+1) * h, y[ix+1]) * (y[ix+1] - y[ix]) - a(t, ix * h, y[ix]) * (y[ix] - y[ix-1])) -
    #       rho(t, ix * h, y[ix]) * y[ix] / tau for ix in range(nx)]

    a = lambda t, ix: 0.5 * (k(t, ix * h, y[ix]) + k(t, (ix-1) * h, y[ix-1]))

    _a = [a(t, ix) * sigma / h**2 for ix in range(nx)]
    _b = [a(t, ix+1) * sigma / h**2 for ix in range(nx)]
    _c = [-((1 / tau) + ((a(t, ix) + a(t, ix+1)) * sigma / h**2)) for ix in range(nx)]
    _f = [-f(t, ix * h, y[ix]) -
          ((1 - sigma) / h**2) * (a(t, ix+1) * (y[ix+1] - y[ix]) - a(t, ix) * (y[ix] - y[ix-1])) -
          1 * y[ix] / tau for ix in range(nx)]

    return _a, _b, _c, _f


def buildMatrix_iter(nx, h, tau, t, sigma, y_s, k, f, y_n):
    a = lambda t, ix: 0.5 * (k(t, ix * h, y_s[ix]) + k(t, (ix-1) * h, y_s[ix-1]))

    _a = [-a(t, ix) / h**2 for ix in range(nx)]
    _b = [-a(t, ix+1) / h**2 for ix in range(nx)]
    _c = [(1 / tau) + ((a(t, ix) + a(t, ix+1)) / h**2) for ix in range(nx)]
    _f = [f(t, ix * h, y_s[ix]) + y_n[ix] / tau for ix in range(nx)]

    return _a, _b, _c, _f

schemes = [buildMatrix_linear, buildMatrix_nonlinear, buildMatrix_iter]


class Calc:
    def __init__(self):
        self.p = {}
        self.it = 0
        self._a = []
        self._b = []
        self._c = []
        self._f = []
        self.y = []
        self.xstep = .0
        self.tstep = .0
        self.nx = 0
        self.sigma = .0
        self.buildMatrix = buildMatrix_nonlinear
        self.y0, self.k, self.f, self.m1, self.m2 = probs[0]
        self.useIter = False

    def x(self, ix):
        return self.xstep * ix

    def t(self, it):
        return self.tstep * it

    def compare(self, y1, y2):
        t = reduce(lambda x, y: x+y, map(lambda x, y: abs(x-y), y1, y2))
        print(t)
        return t > self.p["epsilon"]

    def reset(self):
        self.it = 0
        self.xstep = self.p["xstep"]
        self.tstep = self.p["tstep"]
        self.nx = int(self.p["width"] / self.xstep)

        self.y0, self.k, self.f, self.m1, self.m2 = probs[self.p["prob"]]
        self.buildMatrix = self.p["scheme"]
        self.useIter = self.p["scheme"] == SCH_ITER

        # self.p["epsilon"] = 10e-7 * self.nx
        self.p["epsilon"] = 10e+7 * self.nx

        # При введении правила Рунге по h надо будет пересчитывать с меньшим шагом!
        self.y = [self.y0(ix * self.xstep) for ix in range(self.nx)]

    def calc_sub(self, t, _y, tstep):
        if not self.useIter:
            y = []
            y.extend(_y)
            y.extend([.0, .0])
            _a, _b, _c, _f = self.buildMatrix(self.nx, self.xstep, tstep, t, self.sigma, y, self.k, self.f)
            return TDMA(_a, _b, _c, _f)
        else:
            y = []
            y.extend(_y)
            y_s = []
            y_s.extend(_y)
            for i in range(3):
                y_s.extend([.0, .0])
                _a, _b, _c, _f = buildMatrix_iter(self.nx, self.xstep, tstep, t, self.sigma, y_s, self.k, self.f, y)
                y_s = TDMA(_a, _b, _c, _f)
            return y_s

    def calc(self):
        t = self.t(self.it)

        y1 = self.calc_sub(t, self.y, self.tstep)


        # Runge Rule
        tstep2 = self.tstep * .5

        y2 = self.calc_sub(t, self.y, tstep2)
        t += tstep2
        y2 = self.calc_sub(t, y2, tstep2)


        self.it += 1
        if self.compare(y1, y2):
            self.y = y2
            self.tstep = tstep2
            self.p["tstep"] = self.tstep
            self.it *= 2.0
        else:
            self.y = y1
