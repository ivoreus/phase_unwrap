% this script compiles the library to produce a .mex DLL binary which can be called from Matlab.

% path to include header files
include_path = fullfile(pwd,'include');
include = ['-I' include_path];

% path to source files
src = fullfile(pwd,'src','libsrncp.cpp');

% output path
out_path = fullfile(pwd,'lib');

% build the library
mex('-v', include, src, '-outdir', out_path);

% add the library folder and include folder to matlab path
addpath(out_path);
addpath(include_path);