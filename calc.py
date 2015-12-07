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


class Method:
    def __init__(self):
        self.prec_deg = 1

    def compare(self, y1, y2, prec):
        return max([abs(y1[i] - y2[i]) / (2**self.prec_deg - 1) for i in xrange(len(y1))]) > prec

    def calc_sub(self, nx, xstep, tstep, t, sigma, _y, k, f, m1, m2):
        y = []
        y.extend(_y)
        y.extend([m2(t), m1(t)])
        _a, _b, _c, _f = self.buildMatrix(nx, xstep, tstep, t, sigma, y, k, f)
        return TDMA(_a, _b, _c, _f)


class MethodLinear(Method):
    def __init__(self):
        Method.__init__(self)

    def buildMatrix(self, nx, h, tau, t, sigma, y, k, f):
        _a = [sigma / h**2] * nx
        _b = [sigma / h**2] * nx
        _c = [-((1 / tau) + (2 * sigma / h**2))] * nx
        _f = [-f(t, ix * h, y[ix]) - (1 - sigma) * (y[ix+1] - 2 * y[ix] + y[ix-1]) / h**2 - y[ix] / tau for ix in xrange(nx)]
        return _a, _b, _c, _f


class MethodNonLinear(Method):
    def __init__(self):
        Method.__init__(self)

    def buildMatrix(self, nx, h, tau, t, sigma, y, k, f):
        # _a = [a(t, ix * h, y[ix]) * sigma / h**2 for ix in range(nx)]
        # _b = [a(t, (ix+1) * h, y[ix+1]) * sigma / h**2 for ix in range(nx)]
        # _c = [-((rho(t, ix * h, y[ix]) / tau) + ((a(t, ix * h, y[ix]) + a(t, (ix+1) * h, y[ix+1])) * sigma / h**2)) for ix in range(nx)]
        # _f = [-f(t, ix * h, y[ix]) -
        #       ((1 - sigma) / h**2) * (a(t, (ix+1) * h, y[ix+1]) * (y[ix+1] - y[ix]) - a(t, ix * h, y[ix]) * (y[ix] - y[ix-1])) -
        #       rho(t, ix * h, y[ix]) * y[ix] / tau for ix in range(nx)]

        a = lambda t, ix: 0.5 * (k(t, ix * h, y[ix]) + k(t, (ix-1) * h, y[ix-1]))

        _a = [a(t, ix) * sigma / h**2 for ix in xrange(nx)]
        _b = [a(t, ix+1) * sigma / h**2 for ix in xrange(nx)]
        _c = [-((1 / tau) + ((a(t, ix) + a(t, ix+1)) * sigma / h**2)) for ix in xrange(nx)]
        _f = [-f(t, ix * h, y[ix]) -
              ((1 - sigma) / h**2) * (a(t, ix+1) * (y[ix+1] - y[ix]) - a(t, ix) * (y[ix] - y[ix-1])) -
              1 * y[ix] / tau for ix in xrange(nx)]

        return _a, _b, _c, _f


class MethodIter(Method):
    def __init__(self):
        Method.__init__(self)

    def buildMatrix(self, nx, h, tau, t, sigma, y_s, k, f, y_n):
        a = lambda _t, _ix: 0.5 * (k(_t, _ix * h, y_s[_ix]) + k(_t, (_ix-1) * h, y_s[_ix-1]))

        _a = [-a(t, ix) / h**2 for ix in xrange(nx)]
        _b = [-a(t, ix+1) / h**2 for ix in xrange(nx)]
        _c = [(1 / tau) + ((a(t, ix) + a(t, ix+1)) / h**2) for ix in xrange(nx)]
        _f = [f(t, ix * h, y_s[ix]) + y_n[ix] / tau for ix in xrange(nx)]

        return _a, _b, _c, _f

    def calc_sub(self, nx, xstep, tstep, t, sigma, _y, k, f, m1, m2):
        y = []
        y.extend(_y)
        y_s = []
        y_s.extend(_y)
        for i in xrange(3):
            y_s.extend([m2(t), m1(t)])
            _a, _b, _c, _f = self.buildMatrix(nx, xstep, tstep, t, sigma, y_s, k, f, y)
            y_s = TDMA(_a, _b, _c, _f)
        return y_s


methods = [MethodLinear, MethodNonLinear, MethodIter]


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
        self.y0, self.k, self.f, self.m1, self.m2 = probs[0]
        self.method = methods[0]()

    def x(self, ix):
        return self.xstep * ix

    def t(self, it):
        return self.tstep * it

    def reset(self):
        self.it = 0
        self.xstep = self.p["xstep"]
        self.tstep = self.p["tstep"]
        self.nx = int(self.p["width"] / self.xstep)

        self.y0, self.k, self.f, self.m1, self.m2 = probs[self.p["prob"]]
        self.method = methods[self.p["scheme"]]()

        self.y = [self.y0(ix * self.xstep) for ix in xrange(self.nx)]

    def calc(self):
        t = self.t(self.it)

        y1 = self.method.calc_sub(self.nx, self.xstep, self.tstep, t, self.sigma, self.y, self.k, self.f, self.m1, self.m2)

        # Runge Rule
        tstep2 = self.tstep * .5

        y2 = self.method.calc_sub(self.nx, self.xstep, tstep2, t, self.sigma, self.y, self.k, self.f, self.m1, self.m2)

        t += tstep2
        y2 = self.method.calc_sub(self.nx, self.xstep, tstep2, t, self.sigma, y2, self.k, self.f, self.m1, self.m2)

        self.it += 1
        if self.method.compare(y1, y2, self.p["epsilon"]):
            self.y = y2
            self.tstep = tstep2
            self.p["tstep"] = self.tstep
            self.it *= 2.0
        else:
            self.y = y1
