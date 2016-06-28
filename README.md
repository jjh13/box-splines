# Box-Splines
----

## Introduction
This is the companion code repository to the paper "Fast and exact evaluation
of box  splines via the PP-form". This code does the explicit decomposition
of any box  spline (see ``boxspline.sage``) and will also will the binary tree
representation of said spline (see ``codegen.sage``, there's also code in
that file to export the tree to a C function).

## Dependencies
- SageMath 6.8 (or above, the initial code was written in Sage 6.8, but has been
    working on 7.2) on an Ubuntu machine, this can be installed by following the
    instructions found here: https://help.ubuntu.com/community/SAGE

## Examples
Using the code is fairly simple, start sage via the command ``sage`` -- this
starts an interactive session. At the prompt, load the libraries
```python
load('helpers.sage')
load('boxspline.sage')
load('codegen.sage')
```

Then, creating a box spline object is as simple as:
```python
# Create the courant element
bs = BoxSpline([(1,0),(0,1),(1,1)])
```
By default this is shifted to the origin, you can turn this off by setting
the second argument to ```False```

```python
bs = BoxSpline([(1,0),(0,1),(1,1)], False)
```
The third and fourth arguments shift and re-weight the spline (respectively).

```python
bs = BoxSpline([(1,0),(0,1),(1,1)], False, vector([-1,-1]), 2)
bs.stable_eval((0,0)) == 2
```
Note the use of ``stable_eval`` in the above example, this does not use the
fastest evaluation scheme (it uses a binary search not discussed in the paper).
To create a fast tree evaluation scheme, use the ```PolyTree``` object.
```python
pt = PolyTree(bs)
pt.eval((0,0))
```
Once a ```PolyTree``` object has been created, you can dump the resulting C code
to a text object via ``export_code``
```python
c_code = pt.export_code()
print c_code
```

## License
Copyright &copy; 2016 Joshua Horacsek


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
