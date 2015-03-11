% ------------------------------------------------
% Testing BicubicSR.m

n_channels = 1; % use 1 to only keep first channel or 3 to keep all 3 in RGB

path = '/Users/petervarshavsky/Dropbox/NYU/superresolution/data/meklit/extracted_sequences/large-translation/';
% im_file_1 = strcat(path, 'meklit_lr_0.png');
% im_file_2 = strcat(path, 'meklit_lr_1.png');
im_file_1 = strcat(path, 'meklit1.png');
im_file_2 = strcat(path, 'meklit2.png');

img1 = imread(im_file_1, 'png');
img2 = imread(im_file_2, 'png');
img1 = img1(:,:,1:n_channels); % only keep one channel
img2 = img2(:,:,1:n_channels); % only keep one channel
y1 = double(img1(:));
y2 = double(img2(:));

% 
% imshow(img);
% 
SR_factor = 2; % magnification
[H, W] = size(img1); % Height and Width

options.L = 1;
options.method = 'variational';
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
options.res = 2; % resolution enhancement factor
options.debug = 0;

[newsk, Lambdak,yhat] = LKvar(y1, 2, y2, 0, 1, 1, zeros(3), 1, options);

theta = newsk(1);
sy = newsk(2) * options.res;
sx = newsk(3) * options.res;

% imshow(reshape(yhat, [64,64]));

% Creating an affine transform to test results
tform_matrix = [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
            
tform = affine2d(tform_matrix);
out_img2 = imwarp(img2, tform);
out_img2 = imtranslate(out_img2, [sx, sy]);
% figure, imshow(out_img2);

% Plotting transformed image
figure, subplot(1,3,1), imshow(img1)
subplot(1,3,2), imshow(img2)
subplot(1,3,3), imshow(out_img2)
