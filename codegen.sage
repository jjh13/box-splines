#!/usr/bin/env python
"""
codegen.sage

This is the ancillary code that actually builds a binary tree from a box spline
object. There's also some codegen facilities provided here. See the README for
a quick example.

Copyright © 2016 Joshua Horacsek

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

from sympy.utilities.codegen import codegen
import re

class PolyLeaf(object):
    """
    Represents the leaf node of a tree.
    """
    def __init__(self, bs, region):
        """
        Arguments:
        bs -- the parent box spline object
        region -- the region that this leaf represents, if 0, this is a node on
            the exterior spline, otherwise this shoule be an element from
            bs.get_pp_regions().
        """
        self.bs = bs
        if region == 0:
            self.pp = 0
            self.poly = None
        else:
            self.pp = region['polynomial']
            self.poly = region['polyhedron']

    def eval(self, p):
        """
        Evaluate the current node at p
        """
        return self.pp.subs({x:p[i] for i,x in enumerate(self.bs.x_)})

class PolyNode(object):
    """
    Represents the internal nodes of the binary tree.
    """
    def __init__(self, bs, poly, planes, regions, depth = 0):
        """
        Arguments:
        bs -- the parent box spline object
        poly -- the polyhedron that represents the volume that should be split
            by planes.

        planes -- the planes left to split along

        regions -- the regions that still need to be put into leaves -- this
            should be a sub list of bs.get_pp_regions().
        """
        self.bs = bs

        # Base case
        if self.count_planes(planes) == 1:
            left_region = 0
            right_region = 0

            # Take the only plane we have
            n, dlist, slist = planes[0]
            d = dlist[0]

            self.normal = n
            self.d = d

            # Barf if we have too many planes
            if len(regions) > 2:
                raise AssertionError("This has the wrong amount of planes? Internal bug, please report")

            # Sort the two regions
            for r in regions:
                c = r['center']
                decision = vector(c).dot_product(vector(n)) - d

                if decision <= 0:
                    left_region = r
                else:
                    right_region = r

            self.left_node  = PolyLeaf(bs, left_region)
            self.right_node = PolyLeaf(bs, right_region)
            return

        # Choose plane to split on
        good_plane = 0

        n, dlist, slist = planes[good_plane]
        mid = len(dlist)//2


        # Partition the d and shell lists
        dlist_left, slist_left   = dlist[0:mid] , slist[0:mid]
        dlist_right, slist_right = dlist[mid+1:], slist[mid+1:]
        d     = dlist[mid]
        shell = slist[mid]

        # Keep track of the plane that we used
        self.normal = n
        self.d = d

        # Split the polyhedron based on the chosen plane
        poly_right, poly_left = split_polyhedron(poly, n, d)

        # Partition the polynomial regions
        regions_left = []
        regions_right = []

        for r in regions:
            c = r['center']
            decision = vector(c).dot_product(vector(n)) - d

            if decision <= 0:
                regions_left += [r]
            else:
                regions_right += [r]

        planes_left =  planes[:good_plane] + planes[good_plane+1:]
        planes_right = planes[:good_plane] + planes[good_plane+1:]

        if len(dlist_left) > 0:
            planes_left = planes_left + [(n, dlist_left, slist_left)]

        if len(dlist_right) > 0:
            planes_right = planes_right + [(n, dlist_right, slist_right)]

        planes_left = self.prune_planes(poly_left, planes_left)
        planes_right = self.prune_planes(poly_right, planes_right)

        # This plane was on the shell, so that means
        # exactly side of it should be zero
        if shell:
            if len(regions_left) == 0:
                self.left_node  = PolyLeaf(bs, 0)
                self.right_node = PolyNode(bs, poly_right, planes_right, regions_right, depth+1)
            else:
                self.left_node = PolyNode(bs, poly_left, planes_left, regions_left, depth+1)
                self.right_node = PolyLeaf(bs, 0)
        else:
            if len(planes_left) == 0:
                self.left_node  = PolyLeaf(bs, regions_left[0])
            else:
                self.left_node  = PolyNode(bs, poly_left,  planes_left,  regions_left, depth+1)

            if len(planes_right) == 0:
                self.right_node  = PolyLeaf(bs, regions_right[0])
            else:
                self.right_node = PolyNode(bs, poly_right, planes_right, regions_right, depth+1)

    def eval(self, p):
        """
        Follow the path down the tree, and evaluate the node that corresponds to
        ``p''.
        """

        decision = vector(p).dot_product(vector(self.normal)) - self.d
        if decision <= 0:
            return self.left_node.eval(p)
        return self.right_node.eval(p)

    def prune_planes(self, poly, planes):
        """
        Remove planes from ``planes'' that don't touch poly.
        """
        ret_planes = []

        if poly is None:
            return ret_planes

        for n, dlist, slist in planes:
            new_dlist = []
            new_slist = []
            for (d,s) in zip(dlist,slist):
                if plane_hits_polyhedron(poly, n, d, bs.s_):
                    new_dlist += [d]
                    new_slist += [s]
            if len(new_dlist) > 0:
                ret_planes += [(n, new_dlist, new_slist)]

        return ret_planes

    def count_planes(self, planes):
        """
        Count the planes in the list of grouped planes ``planes''.
        """
        count = 0
        for n, dlist, slist in planes:
            count += len(dlist)
        return count

class PolyTree(object):
    """
    Container for PolyNodes, also provides codegen functionalites.
    """
    def __init__(self, bs):
        """
        Arguments:
            bs -- a box spline object
        """

        self.bs = bs
        self.name = ""

        # Get planes, ROEs and the polyhedron for this bs
        planes = bs._get_grouped_planes_for_eval()
        regions = bs.get_pp_regions()
        poly = bs.get_polyhedron()

        plane_aug = []

        # Augment the planes
        # This attaches whether or not the plane is
        # on the shell of the polyhedron to the plane value
        for n, dlist in planes:
            shell_list = [False]*len(dlist)
            shell_list[0] = shell_list[-1] = True
            plane_aug += [(n, dlist, shell_list)]
        planes = plane_aug

        self.root = PolyNode(bs, poly, planes, regions)

    def eval(self, p):
        """
        Traverse the tree, find the polynomial that correspnds contains ``p'',
        and evaluate it at ``p''.
        """
        return self.root.eval(p)


    def _recurse_code(self, node, depth):
        """
        Recursively generate code.
        """

        # Polynomial leaf, mark the polynomial and emit a call
        if isinstance(node, PolyLeaf):
            if node.pp == 0:
                return '%sreturn 0;' % ("    "*depth)

            region_id = "__pp_r%d__%s__" % (self.id_counter, self.name)
            pp = node.pp
            try:
                pp = pp.full_simplify()
            except:
                pp = node.pp
            self.pieces += [(region_id, pp)]
            self.id_counter += 1

            # Luckily, this retains ordering X_x
            p = "%s(%s)" % (region_id, ', '.join([str(a) for a in pp.variables()]))
            return "%sreturn %s;" % ("    "*depth, p)


        # Otherwise we have a separating node, and we can emit a dot product
        # plane comparison...
        else:
            plane = node.normal
            d = node.d

            lhs = '+'.join([
                            '%s*%s' % (str(v), float(p)) for (v, p) in
                                
                                [(x,y) for x,y in zip(self.bs.x_, plane) if y != 0] ])
            rhs = float(d)
            indent = depth * "    "

            # ... then recurse
            left_code = self._recurse_code(node.left_node, depth + 1)
            right_code = self._recurse_code(node.right_node, depth + 1)

            code = "%sif( %s < %s ) { \n%s \n%s} else { \n%s \n%s}" % (
                indent, lhs, rhs, left_code, indent, right_code, indent
            )
            return code

    def export_code(self, name = "box_spline"):
        """
        Generates C code for the current tree object. Returns a string with
        a C function named ``name'' (by default this is named box_spline).
        """
        self.pieces = []
        self.id_counter = 0

        pieces = []

        tree_code = self._recurse_code(self.root, 1)
        tree_code = "static double %s(%s) {\n%s\n    return 0;\n}\n" %  (
            name,
            ', '.join(['const %s &%s'% ('double', v) for v in self.bs.x_]),
            tree_code
        )
        for p in self.pieces:
            try:
                if len(p[1].variables()) > 0:
                    pieces += [p]
            except:
                pass
        if len(pieces) > 0:
            [(c_name, c_code), (h_name, c_header)] = codegen(pieces, "C", None, header=False, empty=True)
        else:
            c_code = ""

        # Do some popro over the output, we should probably generate this more
        # cleanly....
        c_code = re.sub(r'^#include ".+$', r'', c_code, 0, re.MULTILINE)
        c_code = re.sub(r'^double __pp', r'static inline double __pp', c_code, 0, re.MULTILINE)
        c_code = re.sub(r'double h', r'const double &h', c_code, 0, re.MULTILINE)
        c_code = re.sub(r'double x_([0-9]+)', r'const double &x_\1', c_code, 0, re.MULTILINE)

        return "%s\n%s" % (c_code, tree_code)
