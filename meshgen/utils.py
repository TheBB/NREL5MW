import numpy as np

from math import pi, sin, cos

from GoTools import Point, Surface
import GeoUtils.Interpolate as ip


ex = Point(1, 0, 0)
ey = Point(0, 1, 0)
ez = Point(0, 0, 1)


def add_if_has(obj, attrib, lst):
    if hasattr(obj, attrib):
        lst.append(getattr(obj, attrib))


def surface_curve(srf, k, dirc):


def merge_surfaces(srfa, srfb, dirc):



def mkcircle(center, radius, angle, nelems):
    alpha = angle * pi / 180
    thetas = np.linspace(0, 2*pi, nelems+1)
    pts = [center + Point(radius * cos(t + alpha),
                          radius * sin(t + alpha), 0)
           for t in thetas]

    return ip.CubicCurve(pts=pts, boundary=ip.PERIODIC).ReParametrize(0, 2*pi*radius)


def grading(length, ds, n, tol=1e-12, maxiters=1000):
    fr = lambda s, ds, n, r: s - ds*(1.0 - r**n) / (1.-r)
    dfrdr = lambda ds, n, r: ds * (r**(n-1) * (n*(1.-r) + r) - 1.0) / (1.-r)**2

    # cell size for uniform mesh
    ds_unif = length / n

    # Find start guess for refinement factor
    r = 0.9 if ds_unif <= ds else 1.1

    # Newton loop
    its = 0
    eps = 10 * tol

    for it in xrange(1, maxiters + 1):
        if eps <= tol:
            break

        fval = fr(length, ds, n, r)
        r -= fval / dfrdr(ds, n, r)
        eps  = abs(fval)

    return r if it < maxiters else 0.0


def _F(s, dx1, dx2, n1, n2, r10, r20):
    r2 = r20
    n1i = 1.0 / n1
    coeff1 = dx1 / dx2
    r2n2 = pow(r2, n2)
    r1 = pow(coeff1*r2n2, n1i)
    dx = dx2 / r2n2

    return s - dx * ((1.0-pow(r1,n1+1)) / (1.0-r1) + (1.0-pow(r2,n2+1)) / (1.0-r2))


def _dFdr2(dx1, dx2, n1, n2, r10, r20):
    r2 = r20
    r2n2 = pow(r2, n2)
    r2n2p1 = pow(r2, n2+1)
    c2 = (1.0-r2n2p1) / (1.0-r2)
    dc2dr2 = (-(n2+1) * r2n2 * (1.0-r2) - (1.0-r2n2p1)) / pow(1.0-r2,2)

    n1i = 1.0 / n1
    r1 = pow(dx1 / dx2 * r2n2, n1i)

    r1n1   = pow(r1, n1)
    r1n1p1 = pow(r1, n1+1)
    c1 = (1.0-r1n1p1) / (1.0-r1)
    dc1dr1 = (-(n1+1) * r1n1 * (1.0-r1) - (1.0-r1n1p1)) / pow(1.0-r1,2)
    dr1dr2 = n2 / n1 * pow(dx1/dx2, n1i) * pow(r2, n2/n1-1.0)

    dx = dx2 / r2n2
    dxdr2 = -n2 * dx2 / r2n2p1

    return - dxdr2*(c1 + c2) - dx*(dc1dr1*dr1dr2 + dc2dr2)


def grading_double(length, dx1, dx2, n1, n2, r1=1.1, r2=.9, tol=1e-12, maxiters=200):
    rf = lambda r2: pow(dx1 / dx2 * pow(r2, n2), 1./n1)

    eps = 10.0 * tol
    r1 = rf(r2)

    for it in xrange(1, maxiters + 1):
        if eps <= tol:
            break

        P = _F(length, dx1, dx2, n1, n2, r1, r2)
        dP = _dFdr2(dx1, dx2, n1, n2, r1, r2)
        r2 -= P / dP
        r1 = rf(r2)
        eps = abs(P)

    return (r1, r2) if it < maxiters else (0., 0.)


def gradspace(start, step, factor, N):
    def gen(start, step, factor, N):
        for _ in xrange(N):
            yield start
            start += step
            step *= factor
    return list(gen(float(start), float(step), float(factor), N))
