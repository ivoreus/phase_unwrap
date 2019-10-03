[X,Y,Z] = ndgrid(0:15,0:15,0:15);
im_phase3d_unwrapped = 2*pi*(X*0.2 + Y*0.1 + Z*0.05);
Z_image = 100 * exp(1i * im_phase3d_unwrapped);
snr = 5;
load('randn_save.mat');
%randnRi = randn([16 16 16 16]);
%randnIi = randn([16 16 16 16]);
Ri = real(Z_image); Ri = Ri + (1/snr)*max(abs(Ri(:)))*randnRi(:,:,:,1);
Ii = imag(Z_image); Ii = Ii + (1/snr)*max(abs(Ii(:)))*randnIi(:,:,:,1);
Z_image = Ri + 1i * Ii;
im_phase3d = angle(Z_image);

%% Load library
if not(libisloaded('libsrncp'))
    loadlibrary('libsrncp');
end

im_mask = [];
%% Calculate phase quality map
tic;
[im_unwrapped, im_mask_new] = calllib('libsrncp','Unwrap3D',size(im_phase3d,1),size(im_phase3d,2),size(im_phase3d,3),im_phase3d,im_mask,0);
im_unwrapped = reshape(im_unwrapped,size(im_phase3d));

toc;

unloadlibrary('libsrncp');

figure();imagesc4D(im_phase3d,[-16 16]);colormap(redgreen);
figure();imagesc4D(im_unwrapped,[-16 16]);colormap(redgreen);

