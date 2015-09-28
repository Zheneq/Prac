__author__ = 'Root'

from tdma import *
import math


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
        pass

    def f(self, t, x, u):
        _u = t * math.exp(-(x-5)**2)
        return (math.exp(-(x-5)**2) - (4*x**2 - 40*x + 98))*_u

    def x(self, ix):
        return self.xstep * ix

    def t(self, it):
        return self.tstep * it

    def reset(self):
        self.it = 0
        self.xstep = self.p["xstep"]
        self.tstep = self.p["tstep"]
        self.nx = int(self.p["width"] / self.xstep)

        self.y = [0] * self.nx

    def calc(self):
        self.y.extend([.0, .0])

        t = self.t(self.it)

        self._a = [self.sigma / self.xstep**2] * self.nx
        self._b = [self.sigma / self.xstep**2] * self.nx
        self._c = [(1 / self.tstep) - (2 * self.sigma / self.xstep**2)] * self.nx
        self._f = [-self.f(t, self.x(ix), 0) - (1 - self.sigma) * (self.y[ix+1] - 2 * self.y[ix] + self.y[ix-1]) / self.xstep**2 - self.y[ix] / self.tstep for ix in range(self.nx)]

        self.y = TDMA(self._a, self._b, self._c, self._f)

        self.it += 1