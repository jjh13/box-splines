#!/usr/bin/env python
"""
"""
from sympy.utilities.codegen import codegen

class PolyLeaf(object):
    def __init__(self, bs, region):
        self.bs = bs
        if region == 0:
            self.pp = 0
            self.poly = None
        else:
            self.pp = region['polynomial']
            self.poly = region['polyhedron']

    def eval(self, p):
        return self.pp.subs({x:p[i] for i,x in enumerate(self.bs.x_)})

class PolyNode(object):
    def count_planes(self, planes):
        count = 0
        for n, dlist, slist in planes:
            count += len(dlist)
        return count

    def __init__(self, bs, poly, planes, regions, depth = 0):
        self.bs = bs

        print "We're at depth %d" % depth
        print "Planes: ", planes

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
                print "woah, no, wrong amount of regions!"
                return

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

        print planes
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

        print "We'll be splitting on", n, d

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
            print "We chose to split on the shell"
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
        decision = vector(p).dot_product(vector(self.normal)) - self.d
        if decision <= 0:
            return self.left_node.eval(p)
        return self.right_node.eval(p)

    def prune_planes(self, poly, planes):
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

class PolyTree(object):
    def __init__(self, bs):

        self.bs = bs
        self.name = ""

        # Get planes, ROEs and the polyhedron for this bs
        planes = bs.get_grouped_planes_for_eval()
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
        return self.root.eval(p)


    def recurse_code(self, node, depth, derivative=None):

        # Polynomial leaf, mark the polynomial and emit a call
        if isinstance(node, PolyLeaf):
            if node.pp == 0:
                return '%sreturn 0;' % ("    "*depth)

            if derivative is None:
                region_id = "__pp_r%d__%s__" % (self.id_counter, self.name)
                pp = node.pp
                try:
                    pp = pp.full_simplify()
                except:
                    pp = node.pp
                self.pieces += [(region_id, pp)]
                self.id_counter += 1

                #print region_id
                #print pp
                #print pp.variables()
                # Luckily, this retains ordering X_x
                p = "%s(%s)" % (region_id, ', '.join([str(a) for a in pp.variables()]))
                return "%sreturn %s;" % ("    "*depth, p)
            else:
                dpoly = diff(node.pp, var("x_%d" % derivative))
                try:
                    dpoly = dpoly.full_simplify()
                except:
                    pass
                region_id = "__pp_r%d__%s__d%d__" % (self.id_counter, self.name, derivative)
                self.pieces += [(region_id, dpoly)]
                self.id_counter += 1

                try:
                    # Luckily, this retains ordering X_x
                    if len(dpoly.variables()) > 0:
                        p = "%s(%s)" % (region_id, ', '.join([str(a) for a in dpoly.variables()]))
                        return "%sreturn %s;" % ("    "*depth, p)
                except:
                    pass
                return "%sreturn %s;" % ("    "*depth, dpoly)


        # Otherwise we have a separating node, and we can emit a dot product
        # plane comparison
        else:
            plane = node.normal
            d = node.d

            lhs = '+'.join([
                            '%s*%s' % (str(v), float(p)) for (v, p) in
                                filter(lambda (x,y): y!=0, zip(self.bs.x_, plane))])
            rhs = float(d)
            indent = depth * "    "
            left_code = self.recurse_code(node.left_node, depth + 1, derivative)
            right_code = self.recurse_code(node.right_node, depth + 1, derivative)

            code = "%sif( %s < %s ) { \n%s \n%s} else { \n%s \n%s}" % (
                indent, lhs, rhs, left_code, indent, right_code, indent
            )
            return code

    def export_code(self, sisl = None):
        self.pieces = []
        self.id_counter = 0

        derivative_tree = ""
        pieces = []

        # Generate derivatives if we're generating for SISL
        if sisl is not None:
            for dx_i in xrange(self.bs.s_):
                d_code = self.recurse_code(self.root, 1, dx_i)

                derivative_tree += "%sdouble box_spline_d%d(%s%s) {\n%s\n\treturn 0;\n}\n" %  (
                    "static " if sisl is not None else "",
                    dx_i,
                    "const double &h, " if self.bs.uses_scale() else "",
                    ', '.join(['const %s &%s'% ('double', v) for v in self.bs.x_]),
                    d_code
                )

                for p in self.pieces:
                    try:
                        if len(p[1].variables()) > 0:
                            pieces += [p]
                    except:
                        pass
                self.pieces = []



        tree_code = self.recurse_code(self.root, 1)
        tree_code = "%sdouble box_spline(%s%s) {\n%s\n\treturn 0;\n}\n" %  (
            "static " if sisl is not None else "",
            "const double &h, " if self.bs.uses_scale() else "",
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
            [(c_name, c_code), (h_name, c_header)] = codegen(pieces, "C", None, header=False, empty=False)
        else:
            c_code = ""

        tree_code += derivative_tree

        c_code = re.sub(r'^double __pp', r'static inline double __pp', c_code, 0, re.MULTILINE)
        c_code = re.sub(r'double h', r'const double &h', c_code, 0, re.MULTILINE)
        c_code = re.sub(r'double x_([0-9]+)', r'const double &x_\1', c_code, 0, re.MULTILINE)


        if sisl is not None:
            c_code = re.sub(r'^#include.+$', r'', c_code, 0, re.MULTILINE)
            c_code = "%s\n%s" % (c_code, tree_code)
            c_code = re.sub(r'^', r'        ', c_code, 0, re.MULTILINE)

            deval_1 = ""
            deval_2 = ""

            if self.bs.uses_scale():
                variables =  ', '.join(['(double)p[%s]'% i for i,v in enumerate(self.bs.x_)])
                evaluate_1 = "(double)box_spline(1, %s)" % variables
                evaluate_2 = "(double)box_spline(h, %s)" % variables

                deval_1  = "switch(d){\n"
                deval_1 += ''.join(["                case %d: return (double)box_spline_d%d(1,%s);\n" % (i,i,variables) for i in xrange(self.bs.s_)])
                deval_1 += "                default: break;\n            }"
                deval_2  = "switch(d){\n"
                deval_2 += ''.join(["                case %d: return (double)box_spline_d%d(h,%s);\n" % (i,i,variables) for i in xrange(self.bs.s_)])
                deval_2 += "                default: break;\n            }"

            else:
                variables =  ', '.join(['(double)p[%s]'% i for i,v in enumerate(self.bs.x_)])
                evaluate_1 = "(double)box_spline(%s)" % variables
                evaluate_2 = "(double)box_spline(%s)" % variables
                deval_1  = "switch(d){\n"
                deval_1 += ''.join(["                case %d: return (double)box_spline_d%d(%s);\n" % (i,i,variables) for i in xrange(self.bs.s_)])
                deval_1 += "                default: break;\n            }"
                deval_2  = "switch(d){\n"
                deval_2 += ''.join(["                case %d: return (double)box_spline_d%d(%s);\n" % (i,i,variables) for i in xrange(self.bs.s_)])
                deval_2 += "                default: break;\n            }"



            # Figure out the support for this spline
            p = self.bs.get_polyhedron()
            lattice = lattice_sites_in(p)

            # Convolve the set with itself
            l_prime = set()
            for v in lattice:
                for l in lattice:
                    l_prime.add(tuple(vector(l) + vector(v)))
            lattice = list(l_prime)

            lattice_code = "            Eigen::Matrix<int,%d,%d> mat;\n            mat << \n" % (len(lattice), self.bs.s_)
            for s in lattice:
                lattice_code +="                %s,\n" % (', '.join([str(i) for i in s]))
            lattice_code = lattice_code[:-2] + ";\n"
            lattice_code += "            for(int i = 0; i < mat.rows(); i++) ret.push_back(mat.row(i).transpose());\n "

            format_dict = {
                "spline_name": sisl,
                "spline_xi":str(matrix(self.bs.Xi_).transpose()),
                "spline_name_upper": sisl.upper(),
                "box_tree_code": c_code,
                "phi_code": evaluate_1,
                "dphi_code": deval_1,
                "phi_h_code": evaluate_2,
                "dphi_h_code": deval_2,
                "lattice_code": lattice_code
            }
            return self.get_new_sisl_header().format(**format_dict)
        return "%s\n%s" % (c_code, tree_code)



    def get_new_sisl_header(self):
        return """
/**
 * {spline_name}.hpp
 *
 * Default evaluation code for the spline defined by
{spline_xi}
 * Auto-generated from the Sage worksheet.
 */

#include <vector>
#include <string>
#include <tuple>
#include <sisl/primitives.hpp>
#include <sisl/lattice/base_lattice.hpp>
#include <sisl/basis/basis_function.hpp>

#ifndef __SISL_BOX_{spline_name_upper}__
#define __SISL_BOX_{spline_name_upper}__

namespace sisl {{
    class {spline_name} : public basis_function {{
    private:
{box_tree_code}
    public:
        template <int N>
        static const double phi(const vector &p) {{
            return {phi_code};
        }}
        template <int N>
        static const double phi(const double &h, const vector &p) {{
            return {phi_h_code};
        }}

        template <int N>
        static const double dphi(const vector &p, const int &d) {{
            double ret = 0;
            {dphi_code}
            return 0;
        }}

        template <int N>
        static const double dphi(const double &h, const vector &p, const int &d) {{
            double ret = 0;
            {dphi_h_code}
            return 0;
        }}

        template <int N>
        static std::vector<lattice_site> get_integer_support() {{
            std::vector<lattice_site> ret;
{lattice_code}
            return ret;
        }}

        template<int N, class L, class BF>
        static double convolution_sum(const vector &p, const L *lattice) {{
            /*
             * We have to explicitly call the base method for convolution_sum
             **/
            return basis_function::convolution_sum<N,L,BF>(p, lattice);
        }}

        template<int N, class L, class BF>
        static double convolution_sum_deriv(const vector &p, const L *lattice, const int &component) {{
            return basis_function::convolution_sum_deriv<N,L,BF>(p, lattice, component);
        }}
        template<int N, class L, class BF>
        static vector grad_convolution_sum(const vector &p, const L *lattice) {{
            return basis_function::grad_convolution_sum<N,L,BF>(p, lattice);
        }}
        template<int N, class L, class BF>
        static vector grad_convolution_sum(const vector &p, const L* base, const L **lattices) {{
            return basis_function::grad_convolution_sum<N,L,BF>(p, base, lattices);
        }}
        template<int N, class L, class BF>
        static double convolution_sum_h(const vector &p, const L *lattice, const double &h) {{
            return basis_function::convolution_sum_h<N,L,BF>(p, lattice, h);
        }}
        template<int N, class L, class BF>
        static double convolution_sum_deriv_h(const vector &p, const L *lattice, const int &component, const double &h) {{
            return basis_function::convolution_sum_deriv_h<N,L,BF>(p, lattice, component, h);
        }}
        template<int N, class L, class BF>
        static vector grad_convolution_sum_h(const vector &p, const L *lattice, const double &h) {{
            return basis_function::grad_convolution_sum_h<N,L,BF>(p, lattice, h);
        }}
        template<int N, class L, class BF>
        static vector grad_convolution_sum_h(const vector &p, const L* base, const L **lattices, const double &h) {{
            return basis_function::grad_convolution_sum_h<N,L,BF>(p, base, lattices, h);
        }}
    }};
}}

#endif // __SISL_BOX_{spline_name_upper}__

"""



template = """
/**
 * %s
 * This class defines the abstract base for a spatial filter.
 *
 * @author Joshua Horacsek
 */

#ifndef _SISL2_%s_FILTER_AUTOGEN_H_
#define _SISL2_%s_FILTER_AUTOGEN_H_

#include <sisl2/filter/spatial_filter.h>
#include <tuple>

namespace sisl2 {

template<class T>
class %s_filter : public spatial_filter<T> {
public:
    %s_filter() { }
    ~%s_filter() { }

    virtual std::vector<std::tuple<lattice_site, T>> get_filter() const {
        std::vector<std::tuple<lattice_site, T>> ret;
    	lattice_site l(%d);

        %s
        return ret;
    };
private:
    T m_hScale;
};

}

#endif // _SISL2_%s_FILTER_AUTOGEN_H_
"""

def generate_stencil_code(name, stencil):
    s = len(stencil[0][0])
    h_file = "filter_%s_%d.h" % (name, s)
    bigname = name.upper()
    code = '\n        '.join(["l << %s; ret.push_back(std::make_tuple(l,  %s));" % (",".join([str(st) for st in list(v)]) , float(w)) for (v,w) in stencil])

    code = template % (
        h_file,
        bigname, bigname,
        name, name, name,
        s,
        code,
        bigname
    )
    return (h_file, code)


def generate_filter_code(filter, name, header=True):
    ffilter_template = """
template<class I=double>
class dft_{name}_filter : public frequency_filter<I> {{
public:
    dft_{name}_filter() {{ }}
    ~dft_{name}_filter(){{ }}

    virtual sisl_complex filter(const vector<I> &omega) {{
        return private_filter({freq_variables});
    }}
private:
    {private_filter}
}};
"""

    ffilter_header = """
/**
 * dft_{name}_filter.h
 *
 * Generated from SageMath
 *
 * @author Joshua Horacsek
 **/

#ifndef _SISL2_FOURIER_{bigname}_H_AUTO
#define _SISL2_FOURIER_{bigname}_H_AUTO

#include <Eigen/Dense>
#include <sisl2/primitives.h>
#include <sisl2/frequency/discreet_fourier.h>

namespace sisl2 {{
namespace fourier {{

{all_code}

}}
}}
#endif //_SISL2_FOURIER_{bigname}_H_AUTO
"""
    c_code = codegen(("private_filter", filter), "C", None, header=False, empty=False)[0][1]
    c_code = re.sub(r'^double private_filter', r'static inline double private_filter', c_code, 0, re.MULTILINE)
    c_code = re.sub(r'double h', r'const double &h', c_code, 0, re.MULTILINE)
    c_code = re.sub(r'double x_([0-9]+)', r'const double &x_\1', c_code, 0, re.MULTILINE)

    c_code = re.sub(r'^#include.+$', r'', c_code, 0, re.MULTILINE)
    c_code = re.sub(r'^', r'    ', c_code, 0, re.MULTILINE)

    f = ffilter_template.format(**{
        "name": name,
        "private_filter": c_code,
        "freq_variables": ', '.join([("omega[%d]" % i) for i in xrange(2)])
    })
    if header:
        return ffilter_header.format(**{
            "name":name,
            "bigname":name.upper(),
            "all_code":f
        })
    return f

def generate_filter_code_grouped(gname, filter_list):
    template = """
/**
 * dft_{name}_filter.h
 *
 * Generated from SageMath
 *
 * @author Joshua Horacsek
 **/

#ifndef _SISL2_FOURIER_{bigname}_H_AUTO
#define _SISL2_FOURIER_{bigname}_H_AUTO

#include <Eigen/Dense>
#include <sisl2/primitives.h>
#include <sisl2/frequency/discreet_fourier.h>

namespace sisl2 {{
namespace fourier {{

{all_code}

}}
}}
#endif //_SISL2_FOURIER_{bigname}_H_AUTO
"""
    code = ""
    for Q, name in filter_list:
        code += generate_filter_code(Q, gname+"_"+name, False)

    return template.format(**{
        "name": gname,
        "bigname": gname.upper(),
        "all_code":code
    })
