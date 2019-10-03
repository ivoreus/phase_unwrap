# Multi-dimensional phase unwrapping in C++
This repository contains libraries designed for general purpose multi-dimensional phase unwrapping. They are written in C++ as dynamically loaded libraries to be called from Matlab code.

## SRNCP Phase unwrapping

This library implements the 3D and 4D phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path. 3D algorithm is based on the following publication:

***H.S. Abdul-Rahman, M.A. Gdeisat, D.R. Burton, M.J. Lalor, F. Lilley, and C.J. Moore. Fast and robust three-dimensional best path phase unwrapping algorithm. Applied Optics 2007;46(26):6623-35.***
The 4D algorithm is based on the 3D algorithm.

## Requirements

These libraries and Matlab code examples require Matlab R2017a or later and a compatible C++ compiler (such as from Microsoft Visual Studio 2017). The Matlab MEX C++ compiler is required to be set up prior to compilation.
