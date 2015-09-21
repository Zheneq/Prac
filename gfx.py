#!-*-coding:utf-8-*-
import sys

# import PyQt4 QtCore and QtGui modules
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
import pyqtgraph as pg
import numpy as np
from rk import *
from decimal import Decimal
import math


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

( Ui_MainWindow, QMainWindow ) = uic.loadUiType('window.ui')


class MainWindow(QMainWindow):
    """MainWindow inherits QMainWindow"""

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.gfxBase.addLegend()
        self.ui.gfxPhase1.addLegend()
        self.ui.gfxPhase2.addLegend()
        self.ui.gfxPhase3.addLegend()

        self.p = {}
        self.p["D"] = 0.0
        self.res = []
        self.t = 0.0
        self.y = []

        self.inputs = {}
        self.inputs["resetting"] = [
            (self.ui.txtA, "A"),
            (self.ui.txtB, "B"),
            (self.ui.txtC, "C"),
            (self.ui.txtD, "D"),
            (self.ui.txtQ, "Q"),
            (self.ui.txtEps, "EPS")
        ]
        self.inputs["notresetting"] = [
            (self.ui.txtStep, "step"),
            (self.ui.txtNum, "num")
        ]

        self.reset()

    def __del__(self):
        self.ui = None

    def func1(self, t, yv):
        return 1 - self.p["B"] * yv[0] - yv[0] * yv[1]**2 - self.p["Q"] * yv[0] * yv[1] + yv[2]

    def func2(self, t, yv):
        return self.p["A"] * (yv[0] * yv[1]**2 - yv[1] + self.p["D"])

    def func3(self, t, yv):
        return self.p["C"] * (self.p["Q"] * yv[0] * yv[1] - yv[2])

    def parseinput(self, type):
        inputs = self.inputs[type]

        for t in inputs:
            self.p[t[1]] = t[0].text().toDouble()[0]

        if self.p["step"] == 0.0: self.p["step"] = 0.1

        for t in inputs:
            t[0].setText(str(self.p[t[1]]))

    def reset(self):
        self.y = [0. for i in xrange(3)]
        self.res = [[] for i in xrange(3)]
        self.t = 0.0

        self.parseinput("notresetting")
        self.parseinput("resetting")

        # _p = self.p["B"] - (self.p["D"] + 1)**2 / 3.
        # _q = (-2. * (self.p["D"] + 1)**3 + 9. * (self.p["D"] + 1) * self.p["B"] - 27. * self.p["B"] * self.p["D"])
        # _d = _q**2 - 4/27. * _p**3
        #
        # _w_1 = (.5 * (-_q + math.sqrt(_d)))**(1/3.)
        # _w_2 = (.5 * (-_q - math.sqrt(_d)))**(1/3.)

        self.y[1] = .636523 + self.p["EPS"]

        self.y[0] = (self.p["D"] + 1. - self.y[1])/self.p["B"] + self.p["EPS"]
        self.y[2] = self.p["Q"] * self.y[0] * self.y[1] + self.p["EPS"]

        self.cont()

    def cont(self):
        self.parseinput("notresetting")

        gfx = solve([self.func1, self.func2, self.func3], self.y, self.p["step"], self.p["num"], self.t)
        # map(lambda x, y: x.extend(y), self.res, gfx)
        self.res = gfx

        self.y = map(lambda x: x[-1][1], self.res)
        self.t = self.res[0][-1][0]

        drawable = [np.array(self.res[i]) for i in xrange(len(self.res))]
        self.ui.gfxBase.plotItem.clear()
        self.ui.gfxBase.plotItem.legend.items = []
        for i in xrange(len(drawable)):
            self.ui.gfxBase.plotItem.plot(drawable[i], pen=(i, len(drawable)), name="y"+str(i+1))

        drawable = np.array([[self.res[0][i][1], self.res[1][i][1]] for i in xrange(len(self.res[0]))])
        self.ui.gfxPhase1.plotItem.clear()
        self.ui.gfxPhase1.plotItem.legend.items = []
        self.ui.gfxPhase1.plotItem.plot(drawable, pen=(1,1), name="y2(y1)")

        drawable = np.array([[self.res[1][i][1], self.res[2][i][1]] for i in xrange(len(self.res[0]))])
        self.ui.gfxPhase2.plotItem.clear()
        self.ui.gfxPhase2.plotItem.legend.items = []
        self.ui.gfxPhase2.plotItem.plot(drawable, pen=(1,1), name="y3(y2)")

        drawable = np.array([[self.res[2][i][1], self.res[0][i][1]] for i in xrange(len(self.res[0]))])
        self.ui.gfxPhase3.plotItem.clear()
        self.ui.gfxPhase3.plotItem.legend.items = []
        self.ui.gfxPhase3.plotItem.plot(drawable, pen=(1,1), name="y1(y2)")
# -----------------------------------------------------#


if __name__ == '__main__':
    # create application
    app = QApplication(sys.argv)
    app.setApplicationName('Prac')

    # create widget
    w = MainWindow()
    w.setWindowTitle('Prac')
    w.show()

    # connection
    QObject.connect(app, SIGNAL('lastWindowClosed()'), app, SLOT('quit()'))

    # h = 1e-1
    # numsteps = 1000
    # graph2 = []
    # graph = solve([f_], [1.0], h)

    # numsteps = 250
    # graph1 = solve([f0, f1], [0., 1.], h, numsteps, 0.0)
    # graph2 = [ [[h*i, math.exp(2*h*i) * math.sin(h*i)] for i in xrange(numsteps+1)],
    #            [[h*i, math.exp(h*i) * math.cos(2*h*i)] for i in xrange(numsteps+1)] ]

    # numsteps = 1000
    # graph = solve([g0, g1], [-1., 4.], h, numsteps, 0.0)
    # graph.append([[h*i, (h*i)**2 + 3*(h*i) - 1] for i in xrange(numsteps+1)])
    # graph.append([[h*i, 3*(h*i)**2 + 2*(h*i) + 4] for i in xrange(numsteps+1)])

    # x = [np.random.normal(loc=0., scale=2, size=100)]

    # drawable = [np.array(graph1[i]) for i in xrange(len(graph1))]
    # for i in xrange(len(drawable)):
    #     w.ui.graphicsView.plotItem.plot(drawable[i], pen=None, symbolBrush=(255,0,0), symbolPen='w')
    # drawable = [np.array(graph2[i]) for i in xrange(len(graph2))]
    # for i in xrange(len(drawable)):
    #     w.ui.graphicsView.plotItem.plot(drawable[i], pen=(i, len(drawable)))
    # x.ass = x.sa

    # execute application
    sys.exit(app.exec_())