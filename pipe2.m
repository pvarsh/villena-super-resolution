% Pipeline to construct high resolution image.

[s, mess, messid] = mkdir('tempSR');
handles.LRi = 0;
handles.opt.xtrue = 0;
handles.opt.L = 0;

handles.opt.sx_true = zeros(1,4);
handles.opt.sy_true = zeros(1,4);
handles.opt.theta_true_grad = zeros(1,4);

handles.opt.LR_exe = 0;

h = zeros(9);
h(5,5) = 1;
save('tempSR/h_aux.mat', 'h');

% addpath('functions');

% TODO: add .mat reader. See lines 287 - 355 in SuperResolutionv1.m


%% Load low resolution images
%% TODO: load all in directory
%% TODO: grayscale conversion
num_LR = 4;
for k = 1:num_LR
    im_f_name = strcat('images/lena1_',num2str(k),'.png');
    im = imread(im_f_name);
    if k == 1 % initialize container
        [h,w,chan] = size(im);
        images = zeros(h,w,num_LR);
        images(:,:,k) = im;
    end
    images(:,:,k) = im;
end



handles.opt.L = num_LR;
