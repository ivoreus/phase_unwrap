%% Create a test 4d phase image
[X,Y,Z] = ndgrid(0:15,0:15,0:15);
% im_phase4d_unwrapped matrix has 3 dimensions, and its values obey the
% following rule: im_phase4d_unwrapped(i+1,j+1,k+1,p+1) == 
% 2*pi*(i*0.2 + j*0.1 + k*0.05) for i,j,k = [0..15]
im_phase3d_unwrapped = 2*pi*(X*0.2 + Y*0.1 + Z*0.05);
% Create a complex image with phase as in im_phase4d_unwrapped and
% magnitude of 100
Z_image = 100 * exp(1i * im_phase3d_unwrapped);
% Add some Gaussian noise to the complex image
% Matrices randnRi and randnIi contain normally distributed values
load('randn_save.mat');
%randnRi = randn([16 16 16 16]);
%randnIi = randn([16 16 16 16]);
stdev = 50;
Ri = real(Z_image); Ri = Ri + stdev*randnRi(:,:,:,1);
Ii = imag(Z_image); Ii = Ii + stdev*randnIi(:,:,:,1);
Z_image = Ri + 1i * Ii;
% Calculate the (wrapped) phase of the complex image
im_phase3d = angle(Z_image);

%% Load library
if not(libisloaded('libsrncp'))
    loadlibrary('libsrncp');
end

%% Unwrap phase
% No mask specified: use all phase volume
im_mask = [];
tic;
% On input, provide im_phase3d as 'phase' argument and im_mask as 'mask'
% argument; collect im_unwrapped from 'phase' argument
[im_unwrapped, im_mask_new] = calllib('libsrncp','Unwrap3D',size(im_phase3d,1),size(im_phase3d,2),size(im_phase3d,3),im_phase3d,im_mask,0);
% restore the original shape of im_unwrapped
im_unwrapped = reshape(im_unwrapped,size(im_phase3d));
toc;

unloadlibrary('libsrncp');

%% Show results
% Show original phase
figure();imagesc4D(im_phase4d_unwrapped,[-16 16]);colormap(redgreen);
% Show wrapped phase
figure();imagesc4D(im_phase3d,[-16 16]);colormap(redgreen);
% Show unwrapped phase
figure();imagesc4D(im_unwrapped,[-16 16]);colormap(redgreen);

