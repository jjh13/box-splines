#!/usr/bin/env python
"""
"""

from sage.all import *
from itertools import product, combinations, chain, izip_longest
from operator import mul

load("./helpers.sage")


def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)

class BoxSpline:
    def __init__(self, Xi, centered=True, shift = None, weight=1):

        # Get simple parameters for the box spline
        self.s_ = len(Xi[0])
        self.n_ = len(Xi)

        # Define the fourier and spatial variables for this BS
        self.w_ = [var('w_%d' % i) for i in xrange(self.n_)]
        self.x_ = [var('x_%d' % i) for i in xrange(self.s_)]

        factor = 0 if not centered else 1
        self.c_xi = factor * sum([vector(x) for x in Xi])/2
        self.centered = centered

        if shift:
            self.c_xi -= vector(shift)

        self.Xi_ = Xi[:]
        self.kerXi_ = list(matrix(Xi).kernel().basis())
        self.weight = weight

        # Setup caches for each of these objects
        self.greens_cache = None
        self.differ_cache = None
        self.polytope_cache = None
        self.cached_regions = None
        self.polyhedron_cache = None
        self.polyterm_cache = None
        self.polytermx_cache = None
        self.gt_cache = None

    def warmup(self):
        print "Warming up difference op"
        _ = self.get_differece_operator()
        print "Warming up polyhedron"
        _ = self.get_polyhedron()
        print "Heating the greens function"
        _ = self.decompose_greens()
        print "Grouping planes"
        _ = self.get_grouped_planes_for_eval()
        print "Grouping polyterms"
        _ = self.get_polyterms_w_xform()

    def uses_scale(self):
        return False

    def cleanup(self):
        del self.greens_cache
        del self.differ_cache
        del self.polytope_cache
        del self.cached_regions

    def get_fourier_form(self):
        w = vector([var('w_%d' % i) for i in xrange(self.s_)])
        phi = prod([
                    (1 - exp(-I* w.dot_product(vector(xi))))/(I*w.dot_product(vector(xi)))
                    for xi in self.Xi_
        ])
        ret =  phi * exp(I * self.c_xi.dot_product(w)) * self.weight
        return ret

    def get_centered_fourier(self, use_sinc=False):
        w = vector([var('w_%d' % i) for i in xrange(self.s_)])
        if use_sinc:
            phi = prod([
                w.dot_product(vector(xi))
                for xi in self.Xi_
            ])
        else:
            phi = prod([
                        sin(w.dot_product(vector(xi))/2)/(w.dot_product(vector(xi))/2)
                        for xi in self.Xi_
            ])* self.weight
        return phi
    def get_polyhedron(self):
        """
        Returns a Polyhedron object that represents the support for this
        box-spline.
        """
        if self.polyhedron_cache:
            return self.polyhedron_cache

        diff = self.get_differece_operator()
        self.polyhedron_cache = Polyhedron(vertices=list(set([tuple(v) for v,_ in diff])))
        return self.polyhedron_cache

    def get_differece_operator(self):
        if self.differ_cache:
            return self.differ_cache

        D  = {tuple([0]*self.s_):1}

        for xi in self.Xi_:
            Dp = {}
            for v in D:
                p = tuple(vector(v) + vector(xi))
                Dp[p] = -D[v]

            for v in D:
                if v in Dp:
                    Dp[v] += D[v]
                else:
                    Dp[v] = D[v]
            D = Dp

        self.differ_cache = [(vector(k)-self.c_xi, D[k]) for k in D if D[k] != 0]
        return self.differ_cache



    def simplify_and_split_from_ker(self, v, isolate, gfunc):
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

    def search_nullspace(self, zeros):
        K = matrix(self.kerXi_).transpose()
        N = []
        if sum(zeros) == 0:
            print "Danger case!"
            return vector(self.kerXi_[0])
        for i,_ in filter(lambda (i,z): z != 0, enumerate(zeros)):
            N.append(K[i])
        C = matrix(N).transpose().kernel().basis()[:]
        lc = vector(C[0])
        return matrix(self.kerXi_).transpose() * lc

    def decompose_term(self, gfd):
        (coeff, w) = gfd
        if len(filter(lambda x: x != 0, w)) == self.s_:
            return (gfd)

        # These are the terms that MUST be zero
        constraint = map(lambda x: 1 if x == 0 else 0, w)
        simp = self.search_nullspace(constraint)

        # find any appropriate term to simplify
        sterm = filter(lambda (i,v): v != 0, enumerate(simp))[0][0]

        decomposition = []
        for term in  self.simplify_and_split_from_ker(simp, sterm, gfd):
            decomposition += (self.decompose_term(term))


        return decomposition

    def symbolic_expr(self, gfd):
        (coeff, w) = gfd
        return coeff * prod([1/w_[i]**exp for i,exp in enumerate(w)], 1)

    def decompose_greens(self):
        if self.greens_cache:
            return self.greens_cache

        # The starting point for the simplification is the fourier domain greens function
        nu_alpha = [1] * self.n_


        # If we find repetitions of the vector, change the initial representation
        xi_index = {}
        for i, xi in enumerate(self.Xi_):
            xi = tuple(xi)
            if xi not in xi_index:
                xi_index[xi] = i

        for j, xi in enumerate(self.Xi_):
            xi = tuple(xi)
            if xi in xi_index:
                i = xi_index[xi]
                nu_alpha[i] += 1
                nu_alpha[j] -= 1

        print nu_alpha

        gf = (self.weight, nu_alpha)

        if self.n_ == self.s_:
            return [gf]

        # pick any vector from the kernel to simplify about
        if len([i for i in nu_alpha if i == 0]) == 0:
            # This is the base case
            simp = self.kerXi_[0]
        else:
            constraint = map(lambda x: 1 if x == 0 else 0, nu_alpha)
            simp = self.search_nullspace(constraint)


        # pick a non zero term in that vector
        term = filter(lambda (i,v): v != 0, enumerate(simp))[0][0]

        # recursively decompose terms
        final = []
        for term in  self.simplify_and_split_from_ker(simp, term, gf):
            final += self.decompose_term(term)
        greens = zip(final[::2], final[1::2])

        compacted = {}
        for v,k in greens:
            k = tuple(k)
            compacted[k] = v if k not in compacted else (v + compacted[k])

        self.greens_cache = [(compacted[k], list(k)) for k in compacted]
        return self.greens_cache

    def decompose_greens_2(self):
        gf = (self.weight, [1] * self.n_)

        if self.n_ == self.s_:
            return [gf]

        # pick any vector from the kernel to simplify about
        simp = self.kerXi_[0]

        # pick a non zero term in that vector
        term = filter(lambda (i,v): v != 0, enumerate(simp))[0][0]
        terms = self.simplify_and_split_from_ker(simp, term, gf)

        compacted = {}
        for v, k in terms:
            k = tuple(k)
            compacted[k] = v if k not in compacted else (v + compacted[k])

        new_compacted = {}
        for k in compacted:
            pass

        done = True
        while not done:
            compacted_prime
        return terms

    def calc_knotplanes(self):
        # Calculate the set of knot planes at the origin
        if self.s_ == 1:
            H = list(set(self.Xi_))
        else:
            H = filter(lambda x: len(filter(lambda y: y != 0, x)) > 0, set([ncross_product(nt) for nt in combinations(set(self.Xi_), self.s_ - 1)]))
        H = [vector(v).normalized() for v in H]
        #
        Hprime = []
        Hshell = []

        for plane in H:
            d_list = [0]

            for v in self.Xi_:
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

        Hprime = [(d- vector(self.c_xi)*vector(n), n)  for (d, n) in Hprime]
        Hshell = [(d- vector(self.c_xi)*vector(n), n)  for (d, n) in Hshell]

        return Hprime, Hshell

    def _group_planes(self, planes):
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

    def _recursive_split(self, L, P, depth = 0):
        if depth == 0:
            self.pp_num = 0

        if len(L) == 0:
            if self.pp_num:
                print "hit depth %d: %d" % (depth, self.pp_num)
            self.pp_num += 1
            return [P]

        p, D = L[0]
        mid = len(D)//2

        D_A = D[0:mid]
        D_B = D[mid+1:]
        d = D[mid]

        result = []

        B,A = split_polyhedron(P, p, d)

        # Left
        if A is not None:
            Lnew = L[1:]
            if len(D_A) > 0:
                Lnew += [(p, D_A)]
            result += self._recursive_split(Lnew, A, depth+1)
        # Right
        if B is not None:
            Lnew = L[1:]
            if len(D_B) > 0:
                Lnew += [(p, D_B)]
            result += self._recursive_split(Lnew, B, depth+1)

        return result

    def get_grouped_planes_for_eval(self):
        if self.gt_cache:
            return self.gt_cache

        Hprime, _ = self.calc_knotplanes()
        A = set([(d, tuple(a)) for (d,a) in Hprime])
        B = set([(d, tuple(a)) for (d,a) in _])
        self.gt_cache = self._group_planes(A.union(B))
        return self.gt_cache


    def get_regions(self):
        print "Getting differnce operator"
        differential = self.get_differece_operator()
        print "Got differential"
        Hprime, _ = self.calc_knotplanes()
        print "Got knot planes"

        L = self._group_planes(set([(d, tuple(a)) for (d,a) in Hprime]))
        print "grouped planes"

        poly = Polyhedron(vertices = map(lambda (x,y): list(x), differential))
        print "about to split"
        return self._recursive_split(L, poly)

    def transform_term(self, term):
        (coeffecient, w) = term
        v = [var('v_%d' % i) for i in xrange(self.s_)]
        def make_term(i, k):
            if i == 2:
                return -H(-v[i])*v[i]**(k-1)/factorial(k-1)
            return H(v[i])*v[i]**(k-1)/factorial(k-1)

        itr = filter(lambda (s,m): m !=0, enumerate(w))
        xi_sigma1 = matrix([self.Xi_[i] for i,_ in itr]).transpose().inverse()
        transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())

        subs = xi_sigma1 * vector(self.x_)
        for idx, (sigma, mu) in enumerate(itr):
            transform = transform.substitute(v[idx] == subs[idx])
        return transform

    def poly_term(self, term):
        (coeffecient, w) = term
        v = [var('v_%d' % i) for i in xrange(self.s_)]
        def make_term(i, k):
            return v[i]**(k-1)/factorial(k-1)
        def make_heavy(i):
            return H(v[i])

        itr = filter(lambda (s,m): m !=0, enumerate(w))
        xi_sigma1 = matrix([self.Xi_[i] for i,_ in itr]).transpose().inverse()

        transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())
        heavy = prod([make_heavy(i) for i,(s,m) in enumerate(itr)], 1)

        subs = xi_sigma1 * vector(self.x_)
        for idx, (sigma, mu) in enumerate(itr):
            transform = transform.substitute(v[idx] == subs[idx])
            heavy = heavy.substitute(v[idx] == subs[idx])
        return (transform,heavy)

    def poly_term_w_xi(self, term):
        (coeffecient, w) = term
        v = [var('v_%d' % i) for i in xrange(self.s_)]
        def make_term(i, k):
            return v[i]**(k-1)/factorial(k-1)
        def make_heavy(i):
            return H(v[i])

        itr = filter(lambda (s,m): m !=0, enumerate(w))
        xi_sigma1 = matrix([self.Xi_[i] for i,_ in itr]).transpose().inverse()

        transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())
        heavy = prod([make_heavy(i) for i,(s,m) in enumerate(itr)], 1)

        subs = xi_sigma1 * vector(self.x_)
        for idx, (sigma, mu) in enumerate(itr):
            transform = transform.substitute(v[idx] == subs[idx])
            heavy = heavy.substitute(v[idx] == subs[idx])
        return (transform, heavy, matrix(xi_sigma1, RDF))

    def get_poly(self):
        M = 0
        differential = self.get_differece_operator()
        sgreens = self.decompose_greens()
        sgreens = sum([self.transform_term(t) for t in sgreens])

        for pos, val in differential:
            GT = sgreens
            for i,_ in enumerate(pos):
                GT = GT.substitute(self.x_[i] == self.x_[i] - pos[i])
            M += GT * val
        M = M.expand()
        return M

    def get_polyterms(self):
        if self.polyterm_cache:
            return self.polyterm_cache
        print "decomposing greens"
        greens = self.decompose_greens()
        print "getting polyterms"
        self.polyterm_cache = []
        for (pp,hs) in [self.poly_term(t) for t in greens]:
            self.polyterm_cache += [(pp.full_simplify(), hs)]

        print "ok, got polyterms"
        return self.polyterm_cache

    def get_polyterms_w_xform(self):
        if self.polytermx_cache:
            return self.polytermx_cache
        greens = self.decompose_greens()
        print "getting x polyterms"
        self.polytermx_cache = []
        for (pp,hs,xi) in [self.poly_term_w_xi(t) for t in greens]:
            self.polytermx_cache += [(pp.full_simplify(), hs, xi)]

        print "ok, got x polyterms"
        return self.polytermx_cache

    def get_pp_regions(self):
        if self.cached_regions:
            return self.cached_regions

        # Get auxilary info
        regions = self.get_regions()
        greens = self.decompose_greens()
        differential = self.get_differece_operator()
        polyterms = self.get_polyterms_w_xform()
        # Compute all the polynomials in each region
        summary = []


        stack_xi = polyterms[0][2]
        for (pp, hs, xi) in polyterms[1:]:
            stack_xi = stack_xi.stack(xi)

        for idx, polyhedron in enumerate(regions):
            c = polyhedron.center().n()


            # Compute the pp form of the polynomial
            PP = 0
            print "doing region %d/%d " % (idx, len(regions))
            for jdx, (pos, val) in enumerate(differential):
                posr = vector(pos).n()
                bpp = 0

                pts = stack_xi * (c-vector(pos)).n()
                for idx, p in enumerate(grouper(self.s_, pts)):
                    if any([x <= 0. for x in p]):
                        continue
                    pp, _, _ = polyterms[idx]
                    bpp += pp.subs({self.x_[i]:  self.x_[i] - pos[i] for i,_ in enumerate(pos)})

                PP += val*bpp
            print "Okay!"

            summary.append({
                    'polynomial': PP,
                    'center': c,
                    'polyhedron': polyhedron
                })
        self.cached_regions = summary[:]
        return summary

    def recursive_evaluation(self, pt, X):
        pass

    def stable_eval(self, pt):

        if not self.get_polyhedron().interior_contains(pt):
            return 0

        pt = vector(pt)
        planes = self.get_grouped_planes_for_eval()

        planes_for_polyhedron = []
        planes_to_split = []

        for n, dlist in planes:
            n = vector(n)
            psign = None
            dot_prod = n.dot_product(pt)

            region_idx = [0, len(dlist)]
            left_sign = sign(dot_prod - dlist[0])
            on_plane = False

            # Do a binary search over the dlist
            while region_idx[1] - region_idx[0] > 1:
                new_idx = sum(region_idx)//2
                d = dlist[new_idx]

                if dot_prod - d == 0:  # Add both sides
                    planes_for_polyhedron += [[-dlist[new_idx - 1]] + list(n)]
                    planes_for_polyhedron += [[dlist[new_idx + 1]] + list(-n)]
                    planes_to_split += [(n, d)]
                    on_plane = True
                    break

                if sign(dot_prod - d) == left_sign:
                    region_idx[0] = new_idx
                else:
                    region_idx[1] = new_idx

            if not on_plane:
                planes_for_polyhedron += [[-dlist[region_idx[0]]] + list(n)]
                planes_for_polyhedron += [[ dlist[region_idx[1]]] + list(-n)]

        #build ieqns for polyhedron
        poly = Polyhedron(ieqs=planes_for_polyhedron)

        for n,d in planes_to_split:
            poly, _ = split_polyhedron(poly, n, d)
            if poly is None:
                poly = _

        differential = self.get_differece_operator()

        c = poly.center()
        PP = 0
        for pos, val in differential:
            bpp = 0
            for (pp, hs) in self.get_polyterms():
                hs = hs.subs({self.x_[i]: c[i] - pos[i] for i,_ in enumerate(pos)})
                if hs.n() > 0.5:
                    bpp += pp.subs({self.x_[i]:  pt[i] - pos[i] for i,_ in enumerate(pos)})
            PP += val*bpp
        return PP


    def get_pp_region_for_pt(self, pt):

            if pt not in self.get_polyhedron():
                return 0

            global plist
            global poly
            plist = []
            #print "gg"
            L = self.get_grouped_planes_for_eval()
            #print "gotem"
            pp_num = 0
            poly = None
            planes = {}

            def _is_dead_plane(n,d):
                global plist
                ieqns = [ [-d] + list(n) ] + [[-d] + list(p) for (p,d) in plist]
                return len(Polyhedron(ieqs=ieqns).equations()) > 0

            def _add_plane(poly, n, d):
                if poly is None:
                    return Polyhedron(ieqs=[[-d] + list(n)])
                pprime = Polyhedron(ieqs=list(poly.Hrepresentation()) + [[-d] + list(n)])
                if len(pprime.equations()) > 0:
                    return poly
                return pprime

            def _append_plane(n,d):
                hash = tuple(n)
                nhash = tuple(-vector(n))

                if nhash in planes and -d in planes[nhash]:
                    return
                if hash in planes:
                    if d not in planes[hash]:
                        planes[hash].append(d)
                else:
                    planes[hash] = [d]

            def _recursive_eval(L, x, depth = 0):
                global pp_num
                global plane_list
                global plist
                global poly
                if depth == 0:
                    pp_num = 0

                if len(L) == 0:
                    #print "hit depth %d: %d" % (depth, pp_num)
                    pp_num += 1
                    return

                p, D = L[0]
                mid = len(D)//2

                D_A = D[0:mid]
                D_B = D[mid+1:]
                d = D[mid]

                decision = vector(p).dot_product(x) - d
                if decision <= 0:

                    poly = _add_plane(poly, -vector(p), -d)

                    Lnew = L[1:]
                    if len(D_A) > 0:
                        Lnew += [(p, D_A)]
                    _recursive_eval(Lnew, x, depth+1)
                else:
                    poly = _add_plane(poly, vector(p), d)
                    Lnew = L[1:]
                    if len(D_B) > 0:
                        Lnew += [(p, D_B)]
                    _recursive_eval(Lnew, x, depth+1)
            _recursive_eval(L, vector(pt))

            c = poly.center()
            # if is_degenerate(poly):
            #    print poly, pt, c
            #    print "BARF BARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARFBARF"

            differential = self.get_differece_operator()
            polyterms = self.get_polyterms_w_xform()

            stack_xi = polyterms[0][2]
            for (pp, hs, xi) in polyterms[1:]:
                stack_xi = stack_xi.stack(xi)

            PP = 0
            for pos, val in differential:
                bpp = 0
                pts = stack_xi * (c-vector(pos)).n()

                for idx, p in enumerate(grouper(self.s_, pts)):
                    if any([x <= 0. for x in p]):
                        continue
                    pp, _, _ = polyterms[idx]
                    bpp += pp.subs({self.x_[i]:  self.x_[i] - pos[i] for i,_ in enumerate(pos)})
                PP += val*bpp
            return PP


def autocorr(bs, lattice_f=lambda _: True):
    # Construct the auto correlation BS
    autocorr = BoxSpline(bs.Xi_*2, bs.centered, None, bs.weight)

    # Give us a frequency variable
    w = vector([var("w_%d" % i) for i in xrange(bs.s_)])
    alpha = vector([var("a_%d" % i) for i in xrange(bs.s_)])


    bs_p = bs.get_polyhedron()
    bs_v = matrix([vector(v) for v in bs_p.vertices()]).transpose()

    autocorr_p = autocorr.get_polyhedron()
    autocorr_v = matrix([vector(v) for v in autocorr_p.vertices()]).transpose()

    def nint(v):
        return sign(v)*ceil(abs(v))

    lattice = [x for x in product(*[range(nint(min(list(autocorr_v)[i])), nint(max(list(autocorr_v)[i]))+1 ) for i in xrange(autocorr.s_)]) if ((vector(x) in autocorr_p) and lattice_f(x))]
    stencil = [x for x in product(*[range(nint(min(list(bs_v)[i])), nint(max(list(bs_v)[i]))+1 ) for i in xrange(bs.s_)]) if vector(x) in bs_p]

    auto_corr = [(pt, autocorr.stable_eval(pt)) for pt in lattice]
    return [(vector(pt), v) for (pt,v) in auto_corr if v != 0]
