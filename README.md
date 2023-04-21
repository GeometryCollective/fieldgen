# fieldgen v0.02
--------------------------------------------
**authors:** Keenan Crane, Peter Schröder

![fieldgen](icon.svg)

## About

Given a triangulated surface, `fieldgen` computes the smoothest unit vector
field, or more generally, the smoothest unit _n_-vector field (e.g., _n_=1,2,4
for unit vector, line, and cross fields, respectively).  Singularities,
such as sources and sinks, are automatically placed in locations that allow
the field to achieve optimal smoothness.  Such fields can then be used for a
wide variety of computer graphics and geometry processing tasks such as
surface parameterization, quad meshing, architectural geometry, anisotropic
shading, and texture synthesis.

The code is a reference implementation of the paper

   >Felix Knöppel, Keenan Crane, Ulrich Pinkall, Peter Schröder  
   ["Globally Optimal Direction Fields"](http://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf)  
   SIGGRAPH 2013

This version carefully implements the finite element connection Laplacian as
described in the paper, using Chebyshev expansions to ensure good numerics.  It
also supports some sophisticated features, such as holomorphic/anti-holomorphic
energy, alignment with principal curvature directions, and alignment with the
boundary.

The code itself is somewhat messy research code with a very basic user
interface.  Several other implementations are available, which use a less
sophisticated discretization of the connection Laplacian and do not support
some of the features mentioned above.  However, they may be useful for certain
tasks, or in different build settings.  In particular:

   - The `stripes` code provides a simple version of the algorithm,
     as well as editing of singularities and generation of a field-
     aligned parameterization:

       <http://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/code.zip>

   - There is also an implementation in _Directional_, which is built on
     top of [Eigen](http://eigen.tuxfamily.org), and is hence header-only:

       <https://github.com/avaxman/Directional>

Note that `fieldgen` and stripes are both built on top of the CHOLMOD sparse
direct solver, whereas Directional is built on top of Eigen.  The latter can be
easier to install (since it is header only) but can be significantly slower,
since the Eigen Cholseky solver is much less mature than CHOLMOD.  See also
below for easy install instructions for CHOLMOD.

### Version History

* 0.01 (Sep 1, 2013) — Initial release
* 0.02 (Jun 4, 2019) — Added boundary alignment, OBJ output, command line support

## Installation

fieldgen depends on SuiteSparse, which you can obtain from

   <http://faculty.cse.tamu.edu/davis/suitesparse.html>

On most platforms, SuiteSparse can be installed via standard
package managers.  On Mac OS X / HomeBrew, it can be installed via

   ```brew install homebrew/science/suite-sparse```

To build, you will have to edit the Makefile and set the include/lib
paths accordingly.  Some examples are provided.  Once these paths
have been set, simply type

   ```make```

which (barring any compilation/linker errors) should build an executable
called `fieldviz`.


## Running

Once built, you should be able to run the executable by typing

```./fieldviz data/bunny.obj```

(or specifying a path to any mesh file in OBJ format).  You should
see a window showing the mesh and some information in the upper-left
corner.  Hitting `space` will generate the smoothest field on the surface:

![interface](interface.jpg)

Other commands can be accessed via the keyboard:

| key      | action
| -------- | -----------------------------------------------------------------------------
|  `space` | update field
|    `k/K` | increase/decrease the symmetry degree of the field (1=vector, 2=line, 4=cross)
|    `s/S` | adjust the smoothness energy; -1=holomorphic, 0=Dirichlet, 1=antiholomorphic
|    `t/T` | adjust trade off between smoothness and curvature alignment (if enabled)
|      `c` | toggle curvature alignment
|      `b` | toggle boundary alignment
|      `m` | draw smooth shaded
|      `f` | draw faceted (with wireframe)
|      `*` | show/hide singularities
|      `w` | write solution to `out.obj`
|  `` ` `` | take a screenshot
| `escape` | exit

**Note:** curvature alignment works only when the symmetry degree of the field is 2 or 4.

### Input

`fieldgen` assumes that the input is an [oriented](https://en.wikipedia.org/wiki/Orientability) and [manifold](http://15462.courses.cs.cmu.edu/fall2018/lecture/meshes/slide_013) triangle mesh, with or without boundary.  Meshes should be specified in the [Wavefront OBJ file format](https://en.wikipedia.org/wiki/Wavefront_.obj_file).  Note that the resolution of the input mesh will affect the resolution of the output field, since one vector is computed per vertex.  Also note that extremely poor-quality meshes (e.g., with near-zero angles or triangle areas) might cause problems for field generation, though generally the algorithm is very robust (see in particular Figure 16 of the paper by Knöppel et al, listed above).

### Output

Hitting the `w` key will write the current field (and the mesh) to the file
`out.obj` in the working directory.  Files are written as triangle meshes in
[Wavefront OBJ format](https://en.wikipedia.org/wiki/Wavefront_.obj_file).
They also include a tangent vector field, encoded in comment lines at the end
of the file.  The degree of the field is specified by a line of the form

```degree n```

where (for instance) _n_=1 is a unit vector field, _n_=2 is a
line field, and _n_=4 is a cross field.  Individual vectors
are then specified by lines of the form

```field i x y z```

where `i` is the index of the vertex, and `x` `y` `z` are the three components of the tangent vector.  In the case where these vectors encode an _n_-direction field this vector is just one of the _n_ possible vectors.  The other vectors can be obtained by rotating this one around the corresponding  vertex normal, which is given in the usual `vn` line.  Singularities in the field, which are associated with faces, are indicated by lines

   singularity `i` `s`

where `i` is the index of the triangle, and `s` is the degree of the singularity.  All indices are 1-based rather than 0-based.

## Command Line

A command line version of `fieldgen` is also available, which can be useful when running in batch mode, over a network, or in other situations where OpenGL visualization is not available or appropriate (e.g., bundling `fieldgen` into a plugin).

**To build:** follow the same instructions above, but type

```make commandline```

instead of just `make`.  Doing so should produce an executable called `fieldgen` (rather than `fieldviz`).

**To run:** the only mandatory arguments are paths to the input and output meshes (in OBJ format); by default, `fieldgen` will then compute the smoothest unit vector field.  To get a full set of options, type

```./fieldgen```

which should print the usage string

```usage: ./fieldgen OBJ_INPUT_PATH OBJ_OUTPUT_PATH```  
``` [--degree=n] [--alignToCurvature] [--alignToBoundary] [--s=S] [--t=T]```

The command line options are as follows:

* `degree` — field degree (1=unit vector field, 2=line field, 4=cross field, etc.)
* `alignToCurvature` — align field to principal curvature directions
* `alignToBoundary` — align field to the domain boundary
* `s/S` — controls the smoothness energy; -1=holomorphic, 0=Dirichlet, 1=antiholomorphic
* `t/T` — controls the trade off between smoothness and curvature alignment (if enabled)

Note that enabling boundary alignment will override curvature alignment.

## Source

Much of the source code in this archive is just there to support basic stuff
like loading a mesh, solving a linear system, etc.  The key routines are all in

* `Mesh.cpp`
* `KVecDir.cpp`
* `SectionIntegrals.cpp`

The main routines are


* `Mesh::InitKVecDirData()` — setup
* `Mesh::ComputeSmoothest()` — computes smoothest field
* `Mesh::ComputeSmoothestFixedBoundary()` — computes smoothest field aligned to the boundary
* `Mesh::SmoothestCurvatureAlignment()` — computes curvature-aligned field

An example of how these routines should be called is found in `Viewer::mSmoothField()`.


## License

This code is covered by a standard MIT license.

```
Copyright (c) 2013 Keenan Crane and Peter Schröder.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

