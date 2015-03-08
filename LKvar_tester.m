% ------------------------------------------------
% Testing BicubicSR.m

n_channels = 1; % use 1 to only keep first channel or 3 to keep all 3 in RGB

path = '/Users/petervarshavsky/Dropbox/NYU/superresolution/data/meklit/extracted_sequences/1-px displacement/';
im_file_1 = strcat(path, 'meklit_lr_0.png');
im_file_2 = strcat(path, 'meklit_lr_1.png');
img1 = imread(im_file_1, 'png');
img2 = imread(im_file_2, 'png');
y1 = double(img1(:));
y2 = double(img2(:));
% img = img(:,:,1:n_channels); % only keep one channel
disp(size(y1));
% 
% imshow(img);
% 
%  
SR_factor = 2; % magnification
[H, W] = size(img1); % Height and Width

options.L = 1;
options.res = SR_factor;
options.n = H;
options.m = W;
options.N = options.n;% * SR_factor;
options.M = options.m;% * SR_factor;
options.Real = true;
options.theta = [0,0];
options.theta_init = [0,0];
options.sx_init = [0,0];
options.sx = [0,0];
options.sy_init = [0,0];
options.sy = [0,0];
options.LK_maxit = 35;
options.LK_thr = 1e-4;

[newsk, Lambdak,varargout] = LKvar(y1, 2, y2, 0, 1, 1, zeros(3), 1, options);

% A: 1 1
% H: 1 1
% C: 4096 4096
% x: 4096 1

% options.cancel
% options.ref_cancel