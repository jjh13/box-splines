︠3248b3d8-c7dc-489d-981c-dde19aba5eb3s︠
from sage.all import *
from itertools import product, combinations, chain
from operator import mul
from sage.symbolic.function_factory import function_factory
from itertools import product, combinations, chain, izip_longest

import time
import re
import random
import gc
︡ad90ef43-5af7-4bff-9203-ab0709e3971e︡{"done":true}︡
︠4a4e83ce-252b-4ffa-b9f5-ec4459da778ei︠
%md
## Setup the parameters for this box spline
```use_recursive``` specifies whether or not to use the recursive algorithm to determine the polynomial regions, ```center_box_spline``` 
︡a01aa8a1-37f0-4ffa-a646-ef411b0403c4︡{"done":true,"md":"## Setup the parameters for this box spline\n```use_recursive``` specifies whether or not to use the recursive algorithm to determine the polynomial regions, ```center_box_spline```"}
︠bf5b1b85-b0d4-4c4c-80b9-3c980c762a50s︠
use_recursive = False
center_box_spline = False
fast_polyhedra = True

# Box spline is just a collection of direction vectors
Xi = [
    (1,0),
    (0,1),
    (1,1),
    (2,1)
]

︡8cd1267d-a35d-49f6-9b5c-aedcd462aa9c︡{"done":true}︡
︠6184077b-6cdc-4093-a617-19c599fc3671i︠
%md
Cache some important properties that we'll reuse multiple times later
︡218dba06-6096-415b-a815-ead267390a04︡{"done":true,"md":"Cache some important properties that we'll reuse multiple times later"}
︠0ef753a1-6df6-4422-8813-da1ccd6b8107s︠
Xi_ = matrix(Xi)
s_ = len(Xi[0])
n_ = len(Xi)

box_shift = vector([0]*s_)
if center_box_spline:
    box_shift += sum([vector(_) for _ in Xi])/2

x_ = [var('x_%d' % i) for i in xrange(s_)]
    
# Choose the ring for the polyhedra, AA is exact, RDF is fast but can suffer from numerical issuses
# this shouldn't be an issue, because we don't rely on the polyhedra for each region to be exact,
# we only use their centroid (and really we only need a valid interior point)
PRING = RDF if fast_polyhedra else AA
︡ec5ffc80-f518-4484-a162-1c5684dda5d1︡{"done":true}︡
︠066adf81-4ef6-4bcb-9284-665de34f5df5i︠
%md 
## Spline Geometry
In the next cell, we compute the difference operator for the box spline. We also create a polyhedron with the explicit support for the box spline --- the convex hull of the point distribution for the difference operator gives us the support of the spline.
︡9a1dc264-cf44-4b6e-a4f0-ae44d37f66f8︡{"done":true,"md":"## Spline Geometry\nIn the next cell, we compute the difference operator for the box spline. We also create a polyhedron with the explicit support for the box spline --- the convex hull of the point distribution for the difference operator gives us the support of the spline."}
︠4743b418-8381-4d25-8083-1af280c1bc91s︠

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

def get_differece_operator():
    D  = {tuple([0]*s_):1}
    for xi in Xi:
        Dp = {}
        for v in D:
            p = tuple(vector(v) + vector(xi))
            Dp[p] = -D[v]

        for v in D:
            if v in Dp:
                Dp[v] += D[v]
            else:
                Dp[v] = D[v]
        D = {k: Dp[k] for k in Dp if Dp[k] != 0}

    return [(vector(k)-box_shift, D[k]) for k in D if D[k] != 0]

def calc_knotplanes():
    if s_ == 1:
        H = list(set(Xi))
    else:
        H = filter(lambda x: len(filter(lambda y: y != 0, x)) > 0, set([tuple(ncross_product(nt)) for nt in combinations(set(Xi), s_ - 1)]))

    H = [vector(v).normalized() for v in H]
    Hprime = []
    Hshell = []

    for plane in H:
        d_list = [0]

        for v in Xi_:
            dlp = set(d_list[:])
            d = vector(v).dot_product(vector(plane))
            for dp in dlp:
                d_list.append(d + dp)

        d_list = list(set(d_list))
        d_list.sort()

        min_d, max_d = d_list.pop(0), d_list.pop()
        Hshell += [(min_d, plane), (max_d, plane)]

        for d in d_list:
            Hprime.append((d, plane))

    Hprime = [(d- vector(box_shift)*vector(n), n)  for (d, n) in Hprime]
    Hshell = [(d- vector(box_shift)*vector(n), n)  for (d, n) in Hshell]

    return Hprime, Hshell

def group_planes(planes):
    plane_set = {}
    for d, n in planes:
        tn = tuple(n)
        if tn not in plane_set:
            tn = tuple(-vector(n))

            if tn not in plane_set:
                plane_set[tuple(n)] = set([d])
            else:
                plane_set[tn] = plane_set[tn].union(set([-d]))
        else:
            plane_set[tn] = plane_set[tn].union(set([d]))
    for p in plane_set:
        plane_set[p] = sorted(list(plane_set[p]))
    return  list([(n,plane_set[n]) for n in plane_set])

t0 = time.time()
difference_operator = get_differece_operator()
planes, outside_planes = calc_knotplanes()
Supp = Polyhedron(vertices=list(set([tuple(v) for v,_ in difference_operator])))
t1 = time.time()
geom_time = t1-t0

if s_ == 2 or s_ == 3:
    Supp.plot()
︡adfc922a-fd69-4dd9-ab5b-689b8ed0bcb8︡{"file":{"filename":"/home/user/.sage/temp/project-d80dc9a2-0312-4491-ac51-66924da0b7c6/618/tmp_cpMyHP.svg","show":true,"text":null,"uuid":"6f518673-a3e3-4aa2-84a5-1d9b60f851c1"},"once":false}︡{"done":true}︡
︠882f2e30-ea10-422d-bf05-20e7eac6851fi︠
%md
## Slice geometry
We now slice the support of the box spline by all the knot planes, this produces a list of regions that should each contain a polynomial. This is usually the most costly step --- unless you use the recursive form of evaluation to determine the polynomials within each region.
︡4401cde3-a310-4f05-a84c-5f968621077e︡{"done":true,"md":"## Slice geometry\nWe now slice the support of the box spline by all the knot planes, this produces a list of regions that should each contain a polynomial. This is usually the most costly step --- unless you use the recursive form of evaluation to determine the polynomials within each region."}
︠9d238652-631b-471a-910c-e75846c7ef35s︠
def is_degenerate(P):
    """
    Returns wether or not a polyhedron is degenerate. By degenerate we mean
    not a convex region
    """
    return not P.is_compact() or len(P.equations()) > 0 or P.is_empty()

def split_polyhedron(P, plane, d):
    """
    """
    plane = vector(plane)
    left = [list(x) for x in P.Hrepresentation()] + [[-d] + list(plane)]
    right = [list(x) for x in P.Hrepresentation()] + [[d] + list(-plane)]

    P_l = Polyhedron(ieqs=left, base_ring=PRING)
    if(is_degenerate(P_l)):
        P_l = None

    P_r = Polyhedron(ieqs=right, base_ring=PRING)
    if(is_degenerate(P_r)):
        P_r = None
    return (P_l, P_r)

pp_num = 0
def _recursive_split(L, P, depth = 0):
    global pp_num
    if depth == 0:
            pp_num = 0

    if len(L) == 0:
        if pp_num % 10 == 0:
            print "hit depth %d: %d" % (depth, pp_num)
        pp_num += 1
        return [True, pp_num-1, P]

    p, D = L[0]
    mid = len(D)//2

    # Choose the middle plane, to split along
    D_A = D[0:mid]
    D_B = D[mid+1:]
    d = D[mid]

    result = []

    B,A = split_polyhedron(P, p, d)

    Left = None
    # Left
    if A is not None:
        Lnew = L[1:]
        if len(D_A) > 0:
            Lnew += [(p, D_A)]
        Left = _recursive_split(Lnew, A, depth+1)

    Right = None
    # Right
    if B is not None:
        Lnew = L[1:]
        if len(D_B) > 0:
            Lnew += [(p, D_B)]
        Right = _recursive_split(Lnew, B, depth+1)
    return [False, (p,d), Left, Right]

def reg_tree_to_list(region_tree):
    if region_tree[0] == True: # We're at a leaf
        return [region_tree[2]]

    node_type, plane, Left, Right = region_tree
    ret = []
    # Otherwise we're a node and we need to recurse
    if Left is not None:
        ret += reg_tree_to_list(Left)
    if Right is not None:
        ret += reg_tree_to_list(Right)
    return ret

t0 = time.time()
region_tree = _recursive_split(group_planes(planes), Supp)
regions = reg_tree_to_list(region_tree)
t1 = time.time()

region_time = t1-t0
︡e61ef80b-58a1-4b6c-89d6-2abc2c912ab8︡{"stdout":"hit depth 8: 0"}︡{"stdout":"\nhit depth 8: 10"}︡{"stdout":"\nhit depth 7: 20"}︡{"stdout":"\n"}︡{"done":true}︡
︠03e23b67-7c39-4fee-9509-7edd7b8719c5is︠
%md
## Plot regions
Now that we've decomposed the support into disjoint regions based on the spline's mesh, let's take a look at the structure of the regions.
︡fba368f8-98b8-4b98-afbe-cd8292768b36︡{"done":true,"md":"## Plot regions\nNow that we've decomposed the support into disjoint regions based on the spline's mesh, let's take a look at the structure of the regions."}
︠ef132f7d-e1f6-426b-9209-86636c1570ces︠
if s_ <= 3:
    sum([r.plot() for r in regions])
︡1832075d-6c93-426a-916a-bd24659a19c1︡{"file":{"filename":"/home/user/.sage/temp/project-d80dc9a2-0312-4491-ac51-66924da0b7c6/618/tmp_kCNhng.svg","show":true,"text":null,"uuid":"dd94a6e2-d478-47c8-b8f5-a899de9df9b7"},"once":false}︡{"done":true}︡
︠f602835a-1c81-4610-aed9-c7b88a32b12ci︠
%md
# Create a separable form of the Green's function
This is one of the main contributions of the paper. 
︡9bdfb1ae-3b2f-4daf-b019-5c8215f70d03︡{"done":true,"md":"# Create a separable form of the Green's function\nThis is one of the main contributions of the paper."}
︠ee964738-050c-49a4-8664-487f209c2863s︠

kerXi_ = list(Xi_.kernel().basis())

def simplify_and_split_from_ker(v, isolate, gfunc):
    # A term in the greens function (fourier) can be represented as
    # (coeffecient, [deg(w_1), ..., deg(w_n)])

    # If we simplify with a vector from the kernel
    # v = [v_1, v_2, ..., v_n]
    # choosing the i'th term as the focal point
    # v_i =/= 0
    #
    # Then the simplification proceeds as
    # for each v_k in v such that v_k != 0
    # (-coeffecient*v_k/v_j, [deg(w_1), deg(w_2), ..., deg(w_k)-1, ..., deg(w_i)+1, ..., deg(w_n)])

    (coeffecient, w) = gfunc
    i = isolate
    old_v = len(filter(lambda x: x!=0, w))
    new_terms = []

    for k, v_k in filter(lambda (kk, v_kk): v_kk != 0 and kk != i, enumerate(v)):
        w_prime = w[:]
        w_prime[k] -= 1
        w_prime[i] += 1
        new_terms += [(-coeffecient*v[k]/v[i], w_prime)]

    return new_terms

def search_nullspace(zeros):
    K = matrix(kerXi_).transpose()
    N = []

    # This wouldn't work if sum(zeros) == 0
    # however the proof in the paper tells us this will never happen

    for i,_ in filter(lambda (i,z): z != 0, enumerate(zeros)):
        N.append(K[i])
    C = matrix(N).transpose().kernel().basis()[:]
    lc = vector(C[0])
    return matrix(kerXi_).transpose() * lc

def decompose_term(gfd):
    (coeff, w) = gfd
    if len(filter(lambda x: x != 0, w)) == s_:
        return (gfd)

    # These are the terms that MUST be zero
    constraint = map(lambda x: 1 if x == 0 else 0, w)
    simp = search_nullspace(constraint)

    # find any appropriate term to simplify
    sterm = filter(lambda (i,v): v != 0, enumerate(simp))[0][0]

    decomposition = []
    for term in  simplify_and_split_from_ker(simp, sterm, gfd):
        decomposition += (decompose_term(term))
    return decomposition

def decompose_greens():
    # The starting point for the simplification is the fourier domain greens function
    nu_alpha = [1] * n_

    # If we find repetitions of the vector, change the initial representation
    xi_index = {}
    for i, xi in enumerate(Xi):
        xi = tuple(xi)
        if xi not in xi_index:
            xi_index[xi] = i

    for j, xi in enumerate(Xi_):
        xi = tuple(xi)
        if xi in xi_index:
            i = xi_index[xi]
            nu_alpha[i] += 1
            nu_alpha[j] -= 1

    # This is the initial greens function representation
    gf = (1, nu_alpha)

    # if we have exactly s_ non zero components, we're done
    # and we can just return the greeens function
    if len([_ for _ in nu_alpha if _ != 0]) == s_:
        return [gf]

    # pick the first vector from the ordered kernel to simplify about
    if len([i for i in nu_alpha if i == 0]) == 0:
        # This is the base case
        simp = kerXi_[0]
    else:
        constraint = map(lambda x: 1 if x == 0 else 0, nu_alpha)
        simp = search_nullspace(constraint)

    # pick a non zero term in that vector
    term = filter(lambda (i,v): v != 0, enumerate(simp))[0][0]

    # recursively decompose terms
    final = []
    for term in  simplify_and_split_from_ker(simp, term, gf):
        final += decompose_term(term)
    greens = zip(final[::2], final[1::2])

    # Collect like terms in the final output
    compacted = {}
    for v,k in greens:
        k = tuple(k)
        compacted[k] = v if k not in compacted else (v + compacted[k])

    return [(compacted[k], list(k)) for k in compacted]
t0 = time.time()
greens_function = decompose_greens()
t1 = time.time()
greens_time = t1-t0
︡cc23c271-4ac4-42e9-b96f-8c27f602b26f︡{"done":true}︡
︠94d439a2-b052-4dfa-9bef-7ecff272397ci︠
%md
## Derive Polynomials for Regions
Each region has the a unique polynomial within it, if we know the center of a region $R$ (or any interior point) $\mathbf{c}_R \in R^{\mathrm{o}}$, then we can derive the polynomial within that region using Theorem 4.13
$$
M_\Xi(\mathbf{x}) = \sum_{(b,\mathbf{p})\in S_\Xi} \sum_{(c,\mathbf\alpha) \in P_{n-s}} \frac{b \cdot c}{\left|\det\Xi_{\alpha}\right|} { [\![(\Xi_{\alpha}^{-1}\mathbf{x}) - \mathbf{p}]\!]_{+}}^{(\mathbf{\mu_{{\alpha}}}-\mathbf{1})}
$$
since the term  ${ [\![(\Xi_{\alpha}^{-1}\mathbf{x}) - \mathbf{p}]\!]_{+}}$ is an indicator function, we know that, if $x \in R$ we have
$$
M_\Xi(\mathbf{x}) = P_R(\mathbf{x}) := \sum_{(b,\mathbf{p})\in S_\Xi} \sum_{(c,\mathbf\alpha) \in P_{n-s}} \frac{b \cdot c}{\left|\det\Xi_{\alpha}\right|} { [\![(\Xi_{\alpha}^{-1}\mathbf{c}_R) - \mathbf{p}]\!]_{+}^{\mathbf{0}}} {[\![ \left((\Xi_{\alpha}^{-1}\mathbf{x}) - \mathbf{p} \right) ]\!] }^{(\mathbf{\mu_{{\alpha}}}-\mathbf{1})}. $$
The term ${ [\![(\Xi_{\alpha}^{-1}\mathbf{c}_R) - \mathbf{p}]\!]_{+}}$ acts as an indicator function which informs the sum wether or not a polynomial term from the greens function contributes to the overall polynomial in that region.

Additionally, it is also possible to derive the polynomial for a region via the recursive formulation
$$
(n-s)M_\Xi(x) = \sum_{\xi\in\Xi}t_\xi M_{\Xi\backslash \xi}(\mathbf{x}) + (1 - t_\xi) M_{\Xi\backslash \xi}(\mathbf{x}-\xi)
$$
with $t = \Xi^T(\Xi \Xi^T)^{-1}\mathbf{x}$. Typically, one takes $\mathbf{x}$ as a fixed point in this recursion, but if we allow $x$ to be any $x\in R$ then the following holds
$$
M_\Xi(x) = \frac{1}{n-s}\sum_{\xi\in\Xi}t_\xi M_{\Xi\backslash \xi}(\mathbf{c}_R) + (1 - t_\xi) M_{\Xi\backslash \xi}(\mathbf{c}_R-\xi)
$$
where each $t_\xi$ is dependent on $\mathbf{x}$ via the relationship $t = \Xi^T(\Xi \Xi^T)^{-1}\mathbf{x}$. This gives that
$$
    P_{R}(x) = \frac{1}{n-s}\sum_{\xi\in\Xi}t_\xi M_{\Xi\backslash \xi}(\mathbf{c}_R) + (1 - t_\xi) M_{\Xi\backslash \xi}(\mathbf{c}_R-\xi).
$$
Although, in practice we find this is quite expensive to compute for every $R$, it can be turned on in this notebook by setting ```use_recursive=True```.
︡dd089463-bd55-4035-98d8-198e9b1aba03︡{"done":true,"md":"## Derive Polynomials for Regions\nEach region has the a unique polynomial within it, if we know the center of a region $R$ (or any interior point) $\\mathbf{c}_R \\in R^{\\mathrm{o}}$, then we can derive the polynomial within that region using Theorem 4.13\n$$\nM_\\Xi(\\mathbf{x}) = \\sum_{(b,\\mathbf{p})\\in S_\\Xi} \\sum_{(c,\\mathbf\\alpha) \\in P_{n-s}} \\frac{b \\cdot c}{\\left|\\det\\Xi_{\\alpha}\\right|} { [\\![(\\Xi_{\\alpha}^{-1}\\mathbf{x}) - \\mathbf{p}]\\!]_{+}}^{(\\mathbf{\\mu_{{\\alpha}}}-\\mathbf{1})}\n$$\nsince the term  ${ [\\![(\\Xi_{\\alpha}^{-1}\\mathbf{x}) - \\mathbf{p}]\\!]_{+}}$ is an indicator function, we know that, if $x \\in R$ we have\n$$\nM_\\Xi(\\mathbf{x}) = P_R(\\mathbf{x}) := \\sum_{(b,\\mathbf{p})\\in S_\\Xi} \\sum_{(c,\\mathbf\\alpha) \\in P_{n-s}} \\frac{b \\cdot c}{\\left|\\det\\Xi_{\\alpha}\\right|} { [\\![(\\Xi_{\\alpha}^{-1}\\mathbf{c}_R) - \\mathbf{p}]\\!]_{+}^{\\mathbf{0}}} {[\\![ \\left((\\Xi_{\\alpha}^{-1}\\mathbf{x}) - \\mathbf{p} \\right) ]\\!] }^{(\\mathbf{\\mu_{{\\alpha}}}-\\mathbf{1})}. $$\nThe term ${ [\\![(\\Xi_{\\alpha}^{-1}\\mathbf{c}_R) - \\mathbf{p}]\\!]_{+}}$ acts as an indicator function which informs the sum wether or not a polynomial term from the greens function contributes to the overall polynomial in that region.\n\nAdditionally, it is also possible to derive the polynomial for a region via the recursive formulation\n$$\n(n-s)M_\\Xi(x) = \\sum_{\\xi\\in\\Xi}t_\\xi M_{\\Xi\\backslash \\xi}(\\mathbf{x}) + (1 - t_\\xi) M_{\\Xi\\backslash \\xi}(\\mathbf{x}-\\xi)\n$$\nwith $t = \\Xi^T(\\Xi \\Xi^T)^{-1}\\mathbf{x}$. Typically, one takes $\\mathbf{x}$ as a fixed point in this recursion, but if we allow $x$ to be any $x\\in R$ then the following holds\n$$\nM_\\Xi(x) = \\frac{1}{n-s}\\sum_{\\xi\\in\\Xi}t_\\xi M_{\\Xi\\backslash \\xi}(\\mathbf{c}_R) + (1 - t_\\xi) M_{\\Xi\\backslash \\xi}(\\mathbf{c}_R-\\xi)\n$$\nwhere each $t_\\xi$ is dependent on $\\mathbf{x}$ via the relationship $t = \\Xi^T(\\Xi \\Xi^T)^{-1}\\mathbf{x}$. This gives that\n$$\n    P_{R}(x) = \\frac{1}{n-s}\\sum_{\\xi\\in\\Xi}t_\\xi M_{\\Xi\\backslash \\xi}(\\mathbf{c}_R) + (1 - t_\\xi) M_{\\Xi\\backslash \\xi}(\\mathbf{c}_R-\\xi).\n$$\nAlthough, in practice we find this is quite expensive to compute for every $R$, it can be turned on in this notebook by setting ```use_recursive=True```."}
︠0ff47b23-4000-4e27-9b94-9ccb915ebdc2s︠
def _sgn_fn(self, x, parent=None, algorithm=None):
    return 0 if x < 0 else (1 if x > 0 else 1/2)
H =  function_factory('H', 1, '\\text{H}', evalf_func=_sgn_fn)

def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

def poly_term_w_xi(term):
    (coeffecient, w) = term
    v = [var('v_%d' % i) for i in xrange(s_)]
    def make_term(i, k):
        return v[i]**(k-1)/factorial(k-1)
    def make_heavy(i):
        return H(v[i])

    itr = filter(lambda (s,m): m !=0, enumerate(w))
    xi_sigma1 = matrix(SR, [Xi_[i] for i,_ in itr]).transpose().inverse()

    transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())
    heavy = prod([make_heavy(i) for i,(s,m) in enumerate(itr)], 1)

    subs = xi_sigma1 * vector(x_)
    for idx, (sigma, mu) in enumerate(itr):
        transform = transform.substitute(v[idx] == subs[idx])
        heavy = heavy.substitute(v[idx] == subs[idx])
    return (transform, heavy, matrix(RDF, xi_sigma1))

def get_polyterms_w_xform():
    polytermxs = []
    for (pp,hs,xi) in [poly_term_w_xi(t) for t in greens_function]:
        polytermxs += [(pp.expand(), hs, xi)]
    return polytermxs

def get_pp_regions():
    # Get auxilary info
    greens = greens_function
    differential = difference_operator
    polyterms = get_polyterms_w_xform()

    # Compute all the polynomials in each region
    summary = []

    stack_xi = polyterms[0][2]
    for (pp, hs, xi) in polyterms[1:]:
        stack_xi = stack_xi.stack(xi)

    for idx, polyhedron in enumerate(regions):
        c = polyhedron.center().n()

        # Compute the pp form of the polynomial
        PP = 0
        if idx % 10 == 0:
            print "doing region %d/%d " % (idx, len(regions))
        for jdx, (pos, val) in enumerate(differential):
            posr = vector(pos).n()
            bpp = 0

            pts = stack_xi * (c-vector(pos)).n()
            for idx, p in enumerate(grouper(s_, pts)):
                if any([x <= 0. for x in p]):
                    continue
                pp, _, _ = polyterms[idx]
                try:
                    bpp += pp.subs({x_[i]:  x_[i] - pos[i] for i,_ in enumerate(pos)}).expand()
                except:
                    bpp += pp.subs({x_[i]:  x_[i] - pos[i] for i,_ in enumerate(pos)})

            PP += val*bpp

        summary.append({
                'polynomial': PP,
                'polyhedron': polyhedron
            })
    return summary

def recurse_polynomial(Xi, pt, x_v, dim=2):
    if Xi.is_square():
        if Xi.det() == 0:
            return 0
        if all([ 0<= _ <=1 for _ in Xi.inverse() * (pt + box_shift)]):
            return 1/abs(Xi.det())
        return 0

    Xinv = (Xi*Xi.transpose())
    if Xinv.det() == 0:
        return 0
    Xinv = Xi.transpose()*Xinv.inverse()
    t = Xinv*(x_v + box_shift)
    rv = 0
    for i,xi in enumerate(Xi.columns()):
        Xi_prime = Xi[:].delete_columns([i])
        xi = vector(xi)
        rv += recurse_polynomial(Xi_prime, pt, x_v, dim)*t[i]
        rv += recurse_polynomial(Xi_prime, pt-xi, x_v-xi, dim)*(1-t[i])
    return rv/(len(Xi.columns()) - dim)


def get_pp_regions_recursive():
    summary = []

    for idx, polyhedron in enumerate(regions):
        c = polyhedron.center().n()
        if idx % 10 == 0:
            print "doing region %d/%d " % (idx, len(regions))
        summary.append({
                'polynomial':recurse_polynomial(Xi_.transpose(), c, vector(x_), s_).expand(),
                'polyhedron':polyhedron
        })
    return summary
pp_regions = []

t0 = time.time()
if use_recursive:
    pp_regions = get_pp_regions_recursive()
else:
    pp_regions = get_pp_regions()
t1 = time.time()
pp_time = t1-t0
print pp_time
︡41a2d050-13e4-4257-954d-3c8e214beb5c︡{"stdout":"doing region 0/28 \ndoing region 10/28 \ndoing region 20/28 "}︡{"stdout":"\n"}︡{"stdout":"0.177476167679\n"}︡{"done":true}︡
︠72ce80b3-e7ee-4916-8a1c-9305dba08cebs︠
def find_outer_regions():
    planes = group_planes(outside_planes)
    for i,p in enumerate(pp_regions):
        outside = []
        P = p['polyhedron']
        for n,d_list in planes:
            for d in d_list:
                intersect = Polyhedron(eqns=[[d] + list(n) ], base_ring=PRING).intersection(P)
                if intersect.dim() == s_ - 1:
                    outside += [(n,d)]
        pp_regions[i]['outside'] = outside
    return pp_regions

pp_region_o =  find_outer_regions()
︡630f9c47-44c7-4607-865c-cbc99e58dc98︡{"done":true}︡
︠ae0e482c-dfd2-4332-b6f8-cb5c23bf038bs︠
## 
︡c7ba5eaf-3d4d-43cb-ab59-5d8a6e28da7c︡{"done":true}︡
︠156034b8-cd8d-43bf-ae89-f10b91d9cee1s︠
def fast_eval(x):
    def _fast_eval(x, rt = region_tree):
        if x not in Supp:
            return 0
        if rt[0] == True:
            index = rt[1]
            return pp_regions[index]['polynomial'].subs(
                {x_[i]: x[i] for i in xrange(s_)}
            )
        # We know we're an internal node, so just
        # recurse
        _, (n,d), Left, Right = rt
        if vector(n)*vector(x) - d <= 0:
            return _fast_eval(x, Left)
        return _fast_eval(x, Right)
    if x not in Supp:
        return 0
    return _fast_eval(x)

def recurse_eval(Xi, pt):
    if Xi.is_square():
        if Xi.det() == 0:
            return 0
        if all([ 0<= _ <=1 for _ in Xi.inverse() * (pt + box_shift)]):
            return 1/abs(Xi.det())
        return 0

    Xinv = (Xi*Xi.transpose())
    if Xinv.det() == 0:
        return 0
    Xinv = Xi.transpose()*Xinv.inverse()
    t = Xinv*(pt + box_shift)
    rv = 0
    for i,xi in enumerate(Xi.columns()):
        Xi_prime = Xi[:].delete_columns([i])
        xi = vector(xi)
        rv += recurse_eval(Xi_prime, pt)*t[i]
        rv += recurse_eval(Xi_prime, pt-xi)*(1-t[i])
    return rv/(len(Xi.columns()) - n_)

︡ac5889ea-64b2-401a-8eb0-a76784caa489︡{"done":true}︡
︠6c9aa243-e5ef-4147-96b3-ec6b87d0bd6ds︠
plot3d(lambda x,y: 4*fast_eval((x,y)), (0,4), (0,4))
︡9e246434-1e73-4b9a-85df-fd9599e5afe8︡{"file":{"filename":"dce57b8f-4213-487f-9fe9-0702aed6b909.sage3d","uuid":"dce57b8f-4213-487f-9fe9-0702aed6b909"}}︡{"done":true}︡
︠e855bdf6-7650-44f0-9134-c686495493des︠
print "Total Precompute Time: ", geom_time, region_time, greens_time, pp_time

︡0b00d4e7-fc74-415f-b1fb-07775df89e2c︡{"stdout":"Total Precompute Time:  0.823674917221 0.29637503624 0.000777006149292 0.237161159515\n"}︡{"done":true}︡
︠fd320b59-4ffb-4522-a70f-5c0da1ed3ef4s︠

︡8b127d68-7f8a-4ad0-811a-67bb9acf86a3︡{"done":true}︡
︠0220ac72-2292-4672-a654-ff80db71c1a4s︠

def pascal(n):
    """
    Yield up to row ``n`` of Pascal's triangle, one row at a time.
    The first row is row 0.
    """
    if not isinstance (n, int):
        raise TypeError ('n must be an integer')
    if n < 0:
        raise ValueError ('n must be an integer >= 0')

    for i in xrange(0,n+1):
        yield [bound_ns(i,j) for j in xrange(1,i+1)]

def pascal_latex(n, out=sys.stdout):
    """
    Generate a Pascal triangle for LaTeX embedding.
    Sends output to the file-like object ``out`` (default: sys.stdout).
    """
    out.write('\\begin{tabular}{r%s}\n' % ('c' * (2 * n + 1)))
    for i, row in enumerate(pascal(n)):
        out.write('$n=%d$:& ' % i)
        out.write('   & ' * (n - i))
        out.write(' &    & '.join ('%2d' % coeff for coeff in row))
        out.write('\\\\\\noalign{\\smallskip\\smallskip}\n')
    out.write ('\\end{tabular}\n')

def bound_ns(n,s):
    return 2^n * binomial(n-1,s-1)*Partitions(n,length=s).cardinality()
pascal_latex(int(6))



︡3bd026d6-2baf-4d02-8832-8590774dd3cd︡{"stdout":"\\begin{tabular}{rccccccccccccc}\n$n=0$:&    &    &    &    &    &    & \\\\\\noalign{\\smallskip\\smallskip}\n$n=1$:&    &    &    &    &    &  2\\\\\\noalign{\\smallskip\\smallskip}\n$n=2$:&    &    &    &    &  4 &    &  4\\\\\\noalign{\\smallskip\\smallskip}\n$n=3$:&    &    &    &  8 &    & 16 &    &  8\\\\\\noalign{\\smallskip\\smallskip}\n$n=4$:&    &    & 16 &    & 96 &    & 48 &    & 16\\\\\\noalign{\\smallskip\\smallskip}\n$n=5$:&    & 32 &    & 256 &    & 384 &    & 128 &    & 32\\\\\\noalign{\\smallskip\\smallskip}\n$n=6$:& 64 &    & 960 &    & 1920 &    & 1280 &    & 320 &    & 64\\\\\\noalign{\\smallskip\\smallskip}\n\\end{tabular}\n"}︡{"done":true}
︠6e1eec1e-5c44-4d2d-b1a1-1bcb926bb0a6s︠

bound_ns(3,4)
︡1a439215-5004-40ac-bc6e-7ff47afc6005︡{"stdout":"0\n"}︡{"done":true}︡
︠c5b4f265-e598-471c-b8d2-9e3efb3a3502︠









