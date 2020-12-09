#!/usr/bin/env python
"""
This file contains numerous helpers for the box-spline decomposition code.
Mostly, the helper functions in this file are to do with splitting polyhedral
regions, but there are some hacky combinatorial functions, and some function
definitions in here as well.

Copyright © 2016, 2020 Joshua Horacsek

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
__author__ = "Joshua Horacsek"

from sage.all_cmdline import *
from itertools import product,zip_longest
from sage.symbolic.function_factory import function_factory
import sympy

import six

if six.PY2:
    from itertools import izip_longest
else:
    from itertools import zip_longest

"""
Define some useful symbolic functions, namely a symbolic heaviside and sinc
function.
"""
def _hside(self, x, parent=None, algorithm=None):
    return 0 if x < 0 else (1 if x > 0 else 1/2)
H =  function_factory('H', 1, '\\text{H}', evalf_func=_hside)
sympy.sinc = sympy.Function("sinc")
def _sinc_fn__(self, x, parent=None, algorithm=None):
    return 2*sin(x/2)/x if x != 0 else 1
sinc = function_factory('sinc', 1, '\\text{sinc}', evalf_func=_sinc_fn__,
        conversions= {'mathematica':'sinc', 'sympy':'sinc'})

def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    if six.PY2:
        return izip_longest(*args, fillvalue=fillvalue)
    else:
        return zip_longest(*args, fillvalue=fillvalue)

def plane_hits_polyhedron(poly, n, d, s = 2):
    """
    Checks whether a plane intersects a polyhedron, it also discards
    intersections at a vertex.
    """
    n = vector(n)
    plane = Polyhedron(eqns=[[-d]+list(n)], base_ring=RDF)
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
    
    L = Polyhedron(ieqs=[[-d] + list(plane)], base_ring=RDF)
    R = Polyhedron(ieqs=[[d] + list(-plane)], base_ring=RDF)
#     left = [list(x) for x in P.Hrepresentation()] + [[-d] + list(plane)]
#     right = [list(x) for x in P.Hrepresentation()] + [[d] + list(-plane)]

    try:
        P_l = P.intersection(L)
        if P_l.dim() < P.dim():
            P_l = None
#         P_l = Polyhedron(ieqs=left, base_ring=RDF)
#         if(is_degenerate(P_l)):
#             P_l = None
    except: # This is a lazy way to catch degenerate polytopes
        P_l = None

    try:
        P_r = P.intersection(R)
        if P_r.dim() < P.dim():
            P_r = None
#         P_r = Polyhedron(ieqs=right)
#         if(is_degenerate(P_r)):
#             P_r = None
    except:
        P_r = None
    return (P_l, P_r)

def ncross_product(vectors):
    """
    Generialized cross product. Takes in a list of n-1 vectors in n-dimensions,
    and returns an orthogonal vector.
    """
    vectors = list(vectors)
    dim = len(vectors[0])
    if len(vectors) != dim - 1:
        return None
    rv = [0]*dim
    for i in range(dim):
        v = [0]*dim
        v[i] = 1
        rv[i] = matrix(vectors + [v]).det()
    return tuple(rv)

def is_same_plane(p1, p2, flip=True):
    """
    Checks whether two planes are the same (normalizes and compares the
    direction).
    """
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
