%% Create a test 4d phase image
[X,Y,Z,T] = ndgrid(0:15,0:15,0:15,0:15);
% im_phase4d_unwrapped matrix has 4 dimensions, and its values obey the
% following rule: im_phase4d_unwrapped(i+1,j+1,k+1,p+1) == 
% 2*pi*(i*0.2 + j*0.1 + k*0.05 + p*0.08) for i,j,k,p = [0..15]
im_phase4d_unwrapped = 2*pi*(X*0.2 + Y*0.1 + Z*0.05 + T*0.08);
% Create a complex image with phase as in im_phase4d_unwrapped and
% magnitude of 100
Z_image = 100 * exp(1i * im_phase4d_unwrapped);
% Add some Gaussian noise to the complex image
% Matrices randnRi and randnIi contain normally distributed values
load('randn_save.mat');
%randnRi = randn([16 16 16 16]);
%randnIi = randn([16 16 16 16]);
stdev = 50;
Ri = real(Z_image); Ri = Ri + stdev*randnRi;
Ii = imag(Z_image); Ii = Ii + stdev*randnIi;
Z_image = Ri + 1i * Ii;
% Calculate the (wrapped) phase of the complex image
im_phase4d = angle(Z_image);

%% Load library
if not(libisloaded('libsrncp'))
    loadlibrary('libsrncp');
end

%% Unwrap phase
% No mask specified: use all phase volume
im_mask4d = [];
tic;
% On input, provide im_phase4d as 'phase' argument and im_mask4d as 'mask'
% argument; collect im_unwrapped from 'phase' argument
im_unwrapped = calllib('libsrncp','Unwrap4D',size(im_phase4d,1),size(im_phase4d,2),size(im_phase4d,3),size(im_phase4d,4),im_phase4d,im_mask4d,0);
% restore the original shape of im_unwrapped
im_unwrapped = reshape(im_unwrapped,size(im_phase4d));
toc;

unloadlibrary('libsrncp');

%% Show results
% Show original phase
figure();imagesc4D(im_phase4d_unwrapped,[-16 16]);colormap(redgreen);
% Show wrapped phase
figure();imagesc4D(im_phase4d,[-16 16]);colormap(redgreen);
% Show unwrapped phase
figure();imagesc4D(im_unwrapped,[-16 16]);colormap(redgreen);
