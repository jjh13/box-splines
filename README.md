# Box-Splines
----

## Introduction
This is the companion code repository to the paper "Fast and exact evaluation
of box  splines via the PP-form". This code does the explicit decomposition
of any box  spline (see ``boxspline.sage``) and will also will the binary tree
representation of said spline (see ``codegen.sage``, there's also code in
that file to export the tree to a C function).

There is also a stand-alone [CoCalc](https://cocalc.com/) worksheet included in 
the repository. If you want to test out this work without committing to a full
Sage installation, you can create a free [CoCalc](https://cocalc.com/) account,
then upload the worksheet ```cocalc-worksheet.sagews''' to your account.

## Dependencies
- SageMath 9.2 --- The original code was written for SageMath 6.8, however SageMath
  is under rapid development, and breaks compatability often. Moreover, older versions of
  Sage are difficult to find, and no longer hosted by SageMath. The code has been updated
  to work with Sage 9.2, but hasn't been as thourougly tested as the code for 6.8. If you
  notice any bugs, please open an issue. Sage 9.2 can be installed by following the
  instructions found here: https://doc.sagemath.org/html/en/installation/
  
## Examples
Using the code is fairly simple, start sage via the command ``sage`` -- this
starts an interactive session. At the prompt, load the sage files
```python
load('helpers.sage')
load('boxspline.sage')
load('codegen.sage')
```

Then, creating a box spline object is as simple as choosing the defining
direction vectors:
```python
# Create the courant element
bs = BoxSpline([(1,0),(0,1),(1,1)])
```
By default this is shifted to the origin, you can turn this off by setting
the second argument to ```False```:

```python
bs = BoxSpline([(1,0),(0,1),(1,1)], False)
```
The third and fourth arguments shift and re-weight the spline (respectively):

```python
bs = BoxSpline([(1,0),(0,1),(1,1)], False, vector([-1,-1]), 2)
bs.stable_eval((0,0)) == 1
```
Note the use of ``stable_eval`` in the above example, this does not use the
fastest evaluation scheme (it uses a binary search not discussed in the paper).
To create a fast tree evaluation scheme, use the ```PolyTree``` object.
```python
pt = PolyTree(bs) # This may take a while...
pt.eval((0,0))
```
Once a ```PolyTree``` object has been created, you can dump the resulting C code
to a text object via ``export_code``
```python
c_code = pt.export_code()
print(c_code)
```

I typically start an interactive Jupyter notebook session
``sage --notebook=jupyter``
which opens an interactive notebook style application in your web browser.

A slightly more interesting example is the Hexagonal box spline,
```python
bs = BoxSpline([(1/2,-sqrt(3)/2),(1/2,sqrt(3)/2), (1, 0)])
```

Again, create a binary tree
```python
pt = PolyTree(bs) # This may take a while...
```

If you're in an interactive Jupyter session, plot it with 
```python
plot3d(lambda x,y: pt.eval((x,y)), (-2,2),(-2,2))
```

Again, you can generate code via
```python
print(pt.export_code("hexagonal_1", "C++"))
```
which should give

```C++
#include <math.h>

static inline double __pp_r0____(const double &x_1) {

   double __pp_r0_____result;
   __pp_r0_____result = (2.0/9.0)*sqrt(3)*(2*sqrt(3)*x_1 + 3);
   return __pp_r0_____result;

}

static inline double __pp_r1____(const double &x_0, const double &x_1) {

   double __pp_r1_____result;
   __pp_r1_____result = -2.0/3.0*sqrt(3)*x_0 + (2.0/3.0)*x_1 + (2.0/3.0)*sqrt(3);
   return __pp_r1_____result;

}

static inline double __pp_r2____(const double &x_0, const double &x_1) {

   double __pp_r2_____result;
   __pp_r2_____result = (2.0/9.0)*sqrt(3)*(3*x_0 + sqrt(3)*x_1 + 3);
   return __pp_r2_____result;

}

static inline double __pp_r3____(const double &x_0, const double &x_1) {

   double __pp_r3_____result;
   __pp_r3_____result = -2.0/3.0*sqrt(3)*x_0 - 2.0/3.0*x_1 + (2.0/3.0)*sqrt(3);
   return __pp_r3_____result;

}

static inline double __pp_r4____(const double &x_0, const double &x_1) {

   double __pp_r4_____result;
   __pp_r4_____result = -2.0/9.0*sqrt(3)*(-3*x_0 + sqrt(3)*x_1 - 3);
   return __pp_r4_____result;

}

static inline double __pp_r5____(const double &x_1) {

   double __pp_r5_____result;
   __pp_r5_____result = -2.0/9.0*sqrt(3)*(2*sqrt(3)*x_1 - 3);
   return __pp_r5_____result;

}

static double hexagonal_1(const double &x_0, const double &x_1) {
    if( x_1*1.0 < 0.0 ) { 
        if( x_0*-0.8660254037844386+x_1*0.5 < 0.0 ) { 
            if( x_0*0.8660254037844386+x_1*0.5 < 0.0 ) { 
                if( x_1*1.0 < -0.8660254037844386 ) { 
                    return 0; 
                } else { 
                    return __pp_r0____(x_1); 
                } 
            } else { 
                if( x_0*-0.8660254037844386+x_1*0.5 < -0.8660254037844386 ) { 
                    return 0; 
                } else { 
                    return __pp_r1____(x_0, x_1); 
                } 
            } 
        } else { 
            if( x_0*0.8660254037844386+x_1*0.5 < -0.8660254037844386 ) { 
                return 0; 
            } else { 
                return __pp_r2____(x_0, x_1); 
            } 
        } 
    } else { 
        if( x_0*-0.8660254037844386+x_1*0.5 < 0.8660254037844386 ) { 
            if( x_0*0.8660254037844386+x_1*0.5 < 0.8660254037844386 ) { 
                if( x_1*1.0 < 0.8660254037844386 ) { 
                    if( x_0*-0.8660254037844386+x_1*0.5 < 0.0 ) { 
                        return __pp_r3____(x_0, x_1); 
                    } else { 
                        if( x_0*0.8660254037844386+x_1*0.5 < 0.0 ) { 
                            return __pp_r4____(x_0, x_1); 
                        } else { 
                            return __pp_r5____(x_1); 
                        } 
                    } 
                } else { 
                    return 0; 
                } 
            } else { 
                return 0; 
            } 
        } else { 
            return 0; 
        } 
    }
    return 0;
}
```
This isn't the prettiest code, but it should work, and it should bequite fast on 
modern CPUs.


## License
Copyright &copy; 2016,2020 Joshua Horacsek


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
