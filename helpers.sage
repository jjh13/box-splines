#!/usr/bin/env python
"""
This file contains numerous helpers for the box-spline decomposition code.
Mostly, the helper functions in this file are to do with splitting polyhedral
regions, but there are some hacky combinatorial functions, and some function
definitions in here as well.



"""


from sage.all_cmdline import *
from itertools import product
from sage.symbolic.function_factory import function_factory
import sympy


def _sgn_fn(self, x, parent=None, algorithm=None):
    return 0 if x < 0 else (1 if x > 0 else 1/2)
H =  function_factory('H', 1, '\\text{H}', evalf_func=_sgn_fn)

sympy.sinc = sympy.Function("sinc")
def _sinc_fn__(self, x, parent=None, algorithm=None):
    return 2*sin(x/2)/x if x != 0 else 1
sinc = function_factory('sinc', 1, '\\text{sinc}', evalf_func=_sinc_fn__, conversions= {'mathematica':'sinc', 'sympy':'sinc'})

def plane_hits_polyhedron(poly, n, d, s = 2):
    n = vector(n)
    plane = Polyhedron(eqns=[[-d]+list(n)])
    return plane.intersection(poly).dimension() >= s - 1

def n_choose_rho(n, rho):
    """
    Multivatiate combinations
    return: n!/rho!
    """
    return factorial(n)/prod([factorial(r) for r in rho],1)

def is_degenerate(P):
    """
    Returns wether or not a polyhedron is degenerate.
    """
    return not P.is_compact() or len(P.equations()) > 0 or P.is_empty()

def split_polyhedron(P, plane, d):
    """
    Splits a polyhedron along the plane defined by plane.dot(x,y,z) = d,
    and returns two polyhedra A, B -- one of which may be None.
    """
    plane = vector(plane)
    left = [list(x) for x in P.Hrepresentation()] + [[-d] + list(plane)]
    right = [list(x) for x in P.Hrepresentation()] + [[d] + list(-plane)]

    try:
        P_l = Polyhedron(ieqs=left)
        if(is_degenerate(P_l)):
            P_l = None
    except: # This is a lazy way to catch degenerate polytopes
        P_l = None

    try:
        P_r = Polyhedron(ieqs=right)
        if(is_degenerate(P_r)):
            P_r = None
    except:
        P_r = None
    return (P_l, P_r)

def ncross_product(vectors):
    """
    """
    vectors = list(vectors)
    dim = len(vectors[0])
    if len(vectors) != dim - 1:
        return None
    rv = [0]*dim
    for i in xrange(dim):
        v = [0]*dim
        v[i] = 1
        rv[i] = matrix(vectors + [v]).det()
    return tuple(rv)

def is_same_plane(p1, p2, flip=True):
    p1 = list(p1)
    p2 = list(p2)

    p1[0] = p1[0]*(-1 if flip else 1)

    n1 = ([i for i in p1 if i != 0][0])
    n2 = [i for i in p2 if i != 0][0]

    return all([ i/n1 == j/n2 for i,j in zip(p1, p2)])

def lattice_sites_in(polyhedron, lattice_f = None):
    """
    Returns the lattice sites within a polyhedron
    for a given lattice
    """

    # If there's no lattice
    if lattice_f is None:
        lattice_f = lambda _: true

    nint = lambda v: sign(v)*ceil(abs(v))
    s = polyhedron.dim()
    v = matrix([vector(v) for v in polyhedron.vertices()]).transpose()

    # Get lattice sites that touch the support
    lattice = [x for x in
        product(*[range(nint(min(list(v)[i])), nint(max(list(v)[i]))+1 )
        for i in xrange(s)])
        if vector(x) in polyhedron and lattice_f(x)
    ]
    return lattice

def bcc_lattice_test(v):
    x,y,z = list(v)
    if z % 2:
        if x%2==0 or y%2==0:
            return False
    else:
        if x%2 or y%2:
            return False
    return True

def fcc_lattice_test(v):
    x,y,z = list(v)
    xr = x % 2
    yr = y % 2
    zr = z % 2

    index = xr | (yr  << 1) | (zr << 2)
    lattice_lut = [0, -1, -1, 1, -1, 2, 3, -1]
    return lattice_lut[index] != -1

def filter_stencil(stencil, lattice_f = lambda _: True, rweight=1):
    return [ (s,w*rweight) for (s, w) in stencil if lattice_f(s)]
