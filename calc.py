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
        return 0
        _u = t * math.exp(-(x-5)**2)
        return math.exp(-(x-5)**2) - (4*x**2 - 40*x + 98)*_u

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

        # self.p["epsilon"] = 10e-7 * self.nx
        self.p["epsilon"] = 10e+7 * self.nx

        # self.y = [0] * self.nx
        self.y = [math.exp(-(ix * self.xstep - 5.0)**2) for ix in range(self.nx)]

    def buildMatrix(self, nx, h, tau, t, sigma, y, a, rho, f):
        # _a = [sigma / h**2] * nx
        # _b = [sigma / h**2] * nx
        # _c = [-((1 / tau) + (2 * sigma / h**2))] * nx
        # _f = [-f(t, ix * h, y[ix]) - (1 - sigma) * (y[ix+1] - 2 * y[ix] + y[ix-1]) / h**2 - y[ix] / tau for ix in range(nx)]

        _a = [a(t, ix * h, y[ix]) * sigma / h**2 for ix in range(nx)]
        _b = [a(t, (ix+1) * h, y[ix+1]) * sigma / h**2 for ix in range(nx)]
        _c = [-((rho(t, ix * h, y[ix]) / tau) + ((a(t, ix * h, y[ix]) + a(t, (ix+1) * h, y[ix+1])) * sigma / h**2)) for ix in range(nx)]
        _f = [-f(t, ix * h, y[ix]) -
              ((1 - sigma) / h**2) * (a(t, (ix+1) * h, y[ix+1]) * (y[ix+1] - y[ix]) - a(t, ix * h, y[ix]) * (y[ix] - y[ix-1])) -
              rho(t, ix * h, y[ix]) * y[ix] / tau for ix in range(nx)]

        return _a, _b, _c, _f

    def calc(self):
        self.y.extend([.0, .0])

        t = self.t(self.it)
        self._a, self._b, self._c, self._f = self.buildMatrix(self.nx, self.xstep, self.tstep, t, self.sigma, self.y, lambda t,x,y: 1, lambda t,x,y: 1, self.f)
        y1 = TDMA(self._a, self._b, self._c, self._f)

        # Runge Rule
        tstep2 = self.tstep * .5

        self._a, self._b, self._c, self._f = self.buildMatrix(self.nx, self.xstep, tstep2, t, self.sigma, self.y, lambda t,x,y: 1, lambda t,x,y: 1, self.f)
        y2 = TDMA(self._a, self._b, self._c, self._f)

        y2.extend([.0, .0])
        t += tstep2
        self._a, self._b, self._c, self._f = self.buildMatrix(self.nx, self.xstep, tstep2, t, self.sigma, y2, lambda t,x,y: 1, lambda t,x,y: 1, self.f)

        y2 = TDMA(self._a, self._b, self._c, self._f)


        self.it += 1
        if self.compare(y1, y2):
            self.y = y2
            self.tstep = tstep2
            self.p["tstep"] = self.tstep
            self.it *= 2.0
        else:
            self.y = y1
