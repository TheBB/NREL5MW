import NREL_TFI as nt
import GeoUtils.Interpolate as ip

from pylab import *
from GoTools import *
from numpy import linspace

figure(1)
clf()

idxs = [4,5,7,9,11,14,17]
cols = 'bgrcmykbgrcmykbgrcm'

params = nt.GetParams(nt.default_params)
wingdef = params.wingdef

for sec, c in zip(wingdef, cols):
    if sec.foil == 'cylinder':
        continue
        # cent = sec['chord'] * (0.25 - sec['ao'] + sec['ac'])
        # x = [cent + sec['chord']/2 * sin(t) for t in linspace(0, 2*pi, 50)]
        # y = [sec['chord']/2 * cos(t) for t in linspace(0, 2*pi, 50)]
    else:
        x, y = nt.TrailingEdgeModificationFt(sec.foil, (5.0e-2)/sec.chord)
        x = [(xi-(sec.ao+0.25-sec.ac)) * sec.chord for xi in x]
        y = [yi * sec.chord for yi in y]
        x, y = nt.RotateCoordinates(x, y, Point(0,0,0), sec.theta)
 
        tau1 = Point(x[1]  - x[0], y[1] - y[0], 0).Normalize()
        tau2 = Point(x[-1] - x[-2], y[-1] - y[-2], 0).Normalize()

        pts = [Point(xi,yi,sec.z) for xi,yi in zip(x,y)]
        ws = ip.CubicCurve(pts=pts, boundary=ip.TANGENT, der=[tau1, tau2])
        kns = ws.GetKnots()

        ptste, wingte = nt.TrailingEdgeCurve(ws)
        imid = (len(ptste) - 1) / 2

        s1 = kns[(len(kns) - 1) / 2]
        s2 = kns[-1] - s1
        ds = abs(ptste[1] - ptste[0])
        r1 = nt.ComputeFactor(0.5 * s1, ds, 100)
        r2 = nt.ComputeFactor(0.5 * s2, ds, 100)

        ss = nt.GradedSpace(0, ds, r1, 101)
        ss += nt.GradedSpace(ss[-1], ss[-1]-ss[-2], 1.0/r1, 101)[1:-1] + [s1]
        ss += nt.GradedSpace(ss[-1], ds, r2, 101)[1:]
        ss += nt.GradedSpace(ss[-1], ss[-1]-ss[-2], 1.0/r2, 101)[1:-1] + [kns[-1]]

        pts1 = [ws.Evaluate(s) for s in ss]
        pts2 = ptste[imid:-1] + pts1 + ptste[1:imid+1]
        crv = ip.CubicCurve(pts=pts2, boundary=ip.PERIODIC)

        x = [crv.Evaluate(k)[0] for k in linspace(crv.GetKnots()[0], crv.GetKnots()[-1], 1000)]
        y = [crv.Evaluate(k)[1] for k in linspace(crv.GetKnots()[0], crv.GetKnots()[-1], 1000)]


    x = [-xi for xi in x]
    plot(x, y, color=c)

gca().set_aspect('equal')
savefig('out.pdf', bbox_inches='tight')
show()

