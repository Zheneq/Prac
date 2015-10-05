#!-*-coding:utf-8-*-
import sys

# import PyQt4 QtCore and QtGui modules
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
import pyqtgraph as pg
import numpy as np
from rk import *
from calc import *
from decimal import Decimal
import math

( Ui_MainWindow, QMainWindow ) = uic.loadUiType('window.ui')


class MainWindow(QMainWindow):
    """MainWindow inherits QMainWindow"""

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.gfxBase.addLegend()

        self.calc = Calc()

        self.inputs = {}
        self.inputs["resetting"] = [
            (self.ui.txtWidth, "width"),
            # (self.ui.txtB, "B"),
            # (self.ui.txtC, "C"),
            # (self.ui.txtD, "D"),
            # (self.ui.txtQ, "Q"),
            # (self.ui.txtEps, "EPS")
            (self.ui.txtXStep, "xstep"),
            (self.ui.txtTStep, "tstep"),
        ]
        self.inputs["notresetting"] = [
            (self.ui.txtNum, "num")
        ]

        self.reset()

    def __del__(self):
        self.ui = None

    def parseinput(self, type):
        inputs = self.inputs[type]

        for t in inputs:
            self.calc.p[t[1]] = t[0].text().toDouble()[0]

        if self.calc.p["xstep"] == 0.0: self.calc.p["xstep"] = 0.01
        if self.calc.p["tstep"] == 0.0: self.calc.p["tstep"] = 0.01

        for t in inputs:
            t[0].setText(str(self.calc.p[t[1]]))

        self.calc.sigma = self.ui.sdrImplicitness.value() / 100

    def reset(self):
        self.parseinput("resetting")
        self.parseinput("notresetting")

        self.calc.reset()

        self.cont()

    def cont(self):
        self.parseinput("notresetting")

        for i in range(int(self.calc.p["num"])):
            self.calc.calc()

        drawable = np.array(self.calc.y)
        self.ui.gfxBase.plotItem.clear()
        self.ui.gfxBase.plotItem.legend.items = []
        self.ui.gfxBase.plotItem.plot(drawable, pen=(0, 1), name="y")
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