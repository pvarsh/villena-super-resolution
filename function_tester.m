% ------------------------------------------------
% Testing bicubic

n_channels = 1 % use 1 to only keep first channel or 3 to keep all 3 in RGB

im_file = '/Users/petervarshavsky/Dropbox//NYU/superresolution/data/meklit/meklit.png';
img = imread(im_file, 'png');
img = img(:,:,n_channels); % only keep one channel
imshow(img);

disp('Image size: ');
disp(size(img));

SR_factor = 2 % magnification
[H, W] = size(img); % Height and Width

options.L = 1;
options.res = SR_factor;
options.n = H;
options.m = W;
options.N = options.n * SR_factor;
options.M = options.m * SR_factor;
options.Real = true;

x = BicubicSR(img, options);
disp('x image size: ');
disp(size(x));
xx = reshape(x, 512, 512, n_channels);
imshow(xx);