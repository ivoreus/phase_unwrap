# SRNCP Phase unwrapping
This folder contains the library implementing the 3D and 4D phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path (SRNCP).

## References
* For the 3D algorithm:

***H.S. Abdul-Rahman, M.A. Gdeisat, D.R. Burton, M.J. Lalor, F. Lilley, and C.J. Moore. Fast and robust three-dimensional best path phase unwrapping algorithm. Applied Optics 2007;46(26):6623-35.***

* The 4D algorithm is based on the 3D algorithm.

## Requirements
This library and Matlab code for testing require Matlab R2017a or later and a compatible C++ compiler (such as from Microsoft Visual Studio 2017). The Matlab MEX C++ compiler is required to be set up prior to compilation.
* To set up the compiler from Matlab command prompt:
```
mex -setup cpp
```

## Compilation

* To compile:
1. Navigate to the root directory of the repository:
```
cd /path/to/repo/libsrncp/
```
2. Run compilation script:
```
build
```
3. Compiled library file will appear in `/path/to/repo/libsrncp/lib` folder.
4. To test, nagivate to `test` folder and run 3d or 4d testing script:
```
cd /path/to/repo/libsrncp/test/
libsrncp_test3d
libsrncp_test4d
```

## Use
* Unwrap3D and Unwrap4D functions use **phase** and **mask** arguments in the header for both input and output. On the function call, the **phase** variable should contain the pointer to the wrapped phase array, which is unwrapped after the execution, rewriting original (wrapped) values. The **mask** variable is a pointer to an array that contains 1 for valid and 0 for invalid voxels. During the unwrapping process, neighboring voxels (along each of 3 or 4 dimensions) are grouped together. In the end, voxels that belong to contiguous regions in the image appear in the same group, and those belonging to discontiguous regions appear in a different group. There may be a phase jump between these regions that will not be unwrapped. The algorithm has a flag **mask_largest_unwrapped_group** to show the positions of the voxels in the largest contiguous region.

* Unwrap3D function takes 6 arguments (check also comments in libsrncp.cpp code in `src` folder)
  - **h**: height, number of rows, integer
  - **w**: width, number of columns, integer
  - **d**: depth, number of planes, integer
  - **phase**: pointer to phase array (wrapped on input, unwrapped on output), pointer to array of doubles
  - **mask**: pointer to mask array, pointer to array of integers
  - **mask_largest_unwrapped_group**: flag used to rewrite mask with the voxel positions in the largest contiguous set of unwrapped voxels

* Unwrap4D function takes 7 arguments (check also comments in libsrncp.cpp code in `src` folder)
  - **h**: height, number of rows, integer
  - **w**: width, number of columns, integer
  - **d**: depth, number of planes, integer
  - **s**: number of slabs (3d volumes), the length of the 4th dimension, or the spissitude
  - **phase**: pointer to phase array (wrapped on input, unwrapped on output), pointer to array of doubles
  - **mask**: pointer to mask array, pointer to array of integers
  - **mask_largest_unwrapped_group**: flag used to rewrite mask with the voxel positions in the largest contiguous set of unwrapped voxels

* Arrays should be stored in the column-major order, i.e. in 4d case to the ith row, jth column, kth plane and pth slab in array corresponds the linear index `[i + j*h + k*h*w + p*h*w*d]`.
