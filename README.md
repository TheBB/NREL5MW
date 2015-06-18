# NREL mesh generator

This is a mesh generator for wind turbine blades.  It produces structured isogeometry meshes for use
with IFEM or OpenFOAM.

Dependencies:
- GeoModeler.py

## Invocation

    geoModeler nrel.py [parameters...]

Each parameter can be specified in the form `parameter=value`.  The parameter value will be
automatically coerced to the proper type (string, integer, floating point or boolean).
Additionally, boolean parameters can be set to **true** by giving `parameter` (just the parameter
name) and to **false** by giving `noparameter` (the parameter name with `no` in front).

Parameters can be collected in files.  Supported formats are JSON,

    {
        "param1": value1,
        "param2": value2,
    }

and YAML,

    param1: value1
    param2: value2

Invoke a parameter file using the special `paramfile=<filename>` argument.  These can be mixed
freely with other parameters, and in each case later parameters will override earlier ones.

The generator will perform a rudimentary sanity check on the given parameters to ensure that no
obviously incorrect specifications are given.

Along with the mesh, a `parameters.yaml` file will be written that contains all the parameters
necessary to reproduce the mesh.  This file also contains a timestamp and the exact SHA hash of the
commit that produced the mesh.

## Examples

The `examples` folder has some sample input files.  Also check the default parameter values in
`meshgen/params.py`.

The `wingdefs` folder has sample wing definition files, one for the NREL 5 MW reference blade, and
one for the NACA 0015 symmetric airfoil.

All the airfoils are in the `airfoils` folder.

## Parameters

### Input and output parameters

#### wingdef

Defines the wing. This parameter is a list of dicts, and must be specified in
JSON or YAML. See `examples/nrel_wingdef.yaml` for an example.

Each element in the list specifies a section of the blade with the following attributes:
- *z*: The location along the z-axis.
- *theta*: The twist angle.
- *chord*: The chord length.
- *ac*: The aerodynamic center.
- *ao*: The aerodynamic origin.
- *foil*: The airfoil geometry. Either `cylinder` or the name of a file in the `airfoils` folder.

#### out

The name of the *folder* that the mesh will be written to.  Careful!  It will be created if it does
not exist, and cleared if it does.

#### format

The output format.  Valid values are `IFEM` and `OpenFOAM`.  Note that OpenFOAM output requires
`order=2` and `walldistance` set to `false`.

For IFEM format, a file called `full.xinp` will be produced that includes all the other relevant
files.

#### nprocs

The number of processors to optimize for.  This is *not* the same as `nprocs_mg`.

### Mesh generation

#### debug

Setting this to `true` will produce debug output in a folder called `out`, which will be cleared if
it already exists.

#### nprocs_mg

Number of processors to use for parallel speedup during mesh generation.

#### walldistance

Set to `true` to produce walldistance output.  This may take a long time!  Only for use with `format=IFEM`.

#### mesh_mode

There are three supported mesh modes.
- *2d*: Produces a 2D mesh for a single airfoil.  The wing definition file must contain only one
airfoil.
- *semi3d*: Produces a 2D mesh at every knot along the blade for running semi 3D simulations.  Also
outputs a `beam.g2` file containing the actual beam.
- *3d*: Produces a full 3D mesh cut off at the tip of the blade.

#### order

The spline order of the mesh.  The supported values are 2 (linear), 3 (quadratic) and 4 (cubic).
Note that OpenFOAM format requires linear geometry.

### Physical properties

The `Re` parameter (Reynolds number) helps controlling the radial resolution near the body (see
below under radial resolution).  The script will also output fluid properties (dynamic viscosity and
density) in IFEM mode.  For this, the parameters `rho` (density) and `velocity` are provided.  The
dynamic viscosity *mu* satifies

*mu* = *rho* × *velocity* × *length* / *Re*,

where *length* is the characteristic length, chosen as the *smallest* chordlength defined in the
wing definition file.

### Boundary conditions

The resulting mesh has five or seven boundary sets: inflow, outflow, wing as well as left and right
slipwall.  3D meshes also have hub and antihub sets.  The following parameters specify to which set
the interfaces inbetween should belong.

E.g. to set the edge between *inflow* and *antihub* to belong to the inflow, set `in_anithub=in`.

The interfaces are:
- Inflow and slipwalls: `in_slip`
- Inflow and hub: `in_hub`
- Inflow and antihub: `in_antihub`
- Outflow and slipwalls: `out_slip`
- Outflow and hub: `out_hub`
- Outflow and antihub: `out_antihub`
- Slipwalls and hub: `slip_hub`
- Slipwalls and antihub: `slip_antihub`
- Wing and hub: `wing_hub`
- Wing and antihub: `wing_antihub`

### Trailing edge modification

To produce an O-mesh, a rounded trailing edge modification will be made to each non-cylindrical
airfoil.  The size of this modification is given by the `len_te` parameter.

For computational purposes, the cylindrical airfoils also require a “fake” trailing edge.  The size
of this can be set with the `len_te_cyl_fac` parameter, which should be specificed in terms of
multiples of `len_te`.

### Angular resolution

To adjust the angular resolution within the O-mesh, control the following parameters:
- `n_te`: Number of elements for the trailing edge part of the airfoil.
- `n_back`: Number of elements for the back part of the airfoil.
- `n_front`: Number of elements for the front part of the airfoil.

In total, the sum of these give the number of elements for *half* the O-mesh.  Due to the
eight-patch structure, it is required that the sum of these parameters is divisible by four.

### Radial resolution

The O-mesh consist of an airfoil embedded in a circle embedded in a square.  The radius of the
circle is controlled by the `radius` parameter.  The square will always have a radius of twice that
of the circle, so that the full O-mesh is a square with sidelengths 4×`radius`.

Control the number of elements within the circle with `n_circle`, and the number of elements in the
square with `n_square`.

The element size near the airfoil is controlled by the `Re` parameter, which specifies
(approximately) the Reynold's number of the intended flow, and thus indirectly the size of the
boundary layer, and the parameter `n_bndlayer` which specifies how many elements should be within
the boundary layer.  The boundary layer will have meshlines orthogonal to the body.

In case the TFI algorithm has problems with very concave airfoils, try setting `smoothing` to
`true`.

### Lengthwise resolution

There are four different modes of controlling the lengthwise resolution, given by the `length_mode`
parameter.  The options are:
- *extruded*: Use this option to extrude a single airfoil.  The wing definition file must have only
one airfoil.  The length of the extrusion is given by the `length` parameter.  (This parameter has
no effect in any other length mode, where the length of the blade is determined by the
*z*-components of the sections in the wing definition file.)
- *uniform*: Produces uniformly sized elements.
- *double*: Produces two separate regions of geometrically sized elements, each equally large (in
terms of number of elements).  The parameters `d_join` and `d_tip` give the element sizes at the hub
and tip, respectively.
- *triple*: Produces three separate regions of geometrically sized elements: from the hub to the
join (see below) and two more (as with double) between the join and the tip.  The parameters
`d_join` and `d_tip` give the element sizes at the *join* and the tip, respectively.

In each case, `n_length` gives the total number of elements in the *z*-direction.  In triple mode,
`n_base` determines the number of elements between the hub and the join.  In any other mode,
`n_base` has no effect, and will produce a warning if it is positive.

The *join* is an airfoil where the geometry changes rapidly (e.g. the first non-cylindrical
airfoil).  To prevent a self-intersecting geometry, an additional linear interpolation step is
performed near the join.  To do this, set `join_index` to the index of the offending section in the
wing definition file, and set `join_adds` to some positive integer (roughly, the strength of the
correction.)  A warning will be produced if `join_adds` > 0 in extruded length mode or 2D mesh
mode.

### Extensions

It is often desirable to have a mesh larger than the O-mesh, without having to inflate the size of
the O-mesh or restricting yourself to square meshes.  For this, mesh extension is supported *behind*
(in the outflow direction), *ahead* (in the inflow direction) and on the *sides* (in the slipwall
directions).

Set `behind`, `ahead` or `sides` to a positive number to enable extension.  These parameters give
the extended distance in terms of multiples of `radius`.  The number of elements in each direction
is controlled by the `n_behind`, `n_ahead` and `n_sides` parameters.  The element sizes are uniform
in each case.

### Subdivision

To control the subdivision of the mesh into patches, the following parameters may be used:
- `p_inner`: The number of patches radially in the O-mesh.
- `p_behind`: The number of patches in the outflow extension.
- `p_ahead`: The number of patches in the inflow extension.
- `p_sides`: The number of patches in the slipwalls extension.
- `p_length`: The number of patches lengthwise.

For FSI purposes, the the innermost patches in the O-mesh can be added to a topologyset called
*rigid*.  The number of patches can be controlled with `p_rigid`.
