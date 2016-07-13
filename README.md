# image2mesh_cgal
Uses CGAL library to convert a stack of 2D images to tetrahedral meshes.

The `wip` branch contains improved version of the code that uses `json` format to read various input parameters that control the mesh generation process.
In order to use this improved version, MATLAB scripts that read and prepare data should be also changes, WHICH IS NOT DONE YET!

The `release` branch has the code that's currently compiled to various binaries and is shipped with nirfast.

To build the `wip` branch, you will need cmake, CGAL, mpfr and gmp. On OSX and Linux, building should be straightforward,
however Windows is a different story!

### For Linux, OSX
- After updating your system's package manager. Create a `build` folder and then run cmake as usual:

```shell
mkdir build
cd build
cmake ..
```


