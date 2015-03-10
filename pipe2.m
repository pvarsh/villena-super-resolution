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
        handles.opt.m = h;
        handles.opt.n = w;
    end
    images(:,:,k) = im;
end

y = []; % all lr_images as column vectors
for k = 1:num_LR
    yk = double(images(:,:,k));
    % ys = yk/(max(max(yk))); % normalize

    yvecs{k} = yk(:);
    y = [y; yk(:)];

end

handles.yvecs = yvecs;

handles.opt.L = num_LR;
handles.opt.res = 2;

%% Pre-LK registration options

% Load Base Options
%% Specific to Bayesian methods
%opt.noise_estimate = 'COMMON';
% handles.opt.noise_estimate = 'SEPARATE';
handles.opt.LK_maxit = 35;
handles.opt.LK_thr = 1e-4;
% handles.opt.PCG_minit = 20;
% handles.opt.PCG_maxit = 100;
% handles.opt.PCG_thr = 1e-10;
% handles.opt.maxit = 100;%opt.maxit = 100;
% handles.opt.thr = 1e-5;

% handles.opt.InitWithAvgImg = 1;
% handles.opt.InitImgExists = 0;

% %% Output options
% handles.opt.WriteImages = 0;
% handles.opt.ShowImages = 1;
% handles.opt.KeepHistory = 0;
% handles.opt.verbose = 1; 
handles.opt.debug = 0;
% handles.opt.interv_ShowImages = 1; % Show images each opt.interv_ShowImages iterations
% handles.opt.FixedBeta = 0;        
%----------------

% set(handles.text_Runing,'Visible','Off');
% set(handles.text_Runing,'Enable','Off');
% set(handles.opt.ref_cancel,'Value',0);
% set(handles.pushbutton_saveSR,'Visible','off');
% set(handles.pushbutton_saveSR,'Enable','off');

% set(handles.pushbutton_BestResults,'Visible','off');
% set(handles.pushbutton_BestResults,'Enable','off');

handles.opt.M =  handles.opt.m;
handles.opt.N =  handles.opt.n; % For the initial registration 

handles.opt.xtrue = [];
handles.opt.method ='variational';
handles.opt.Real = 1;

handles.opt.theta_init = zeros(1, handles.opt.L);
handles.opt.sx_init = zeros(1, handles.opt.L);
handles.opt.sy_init = zeros(1, handles.opt.L);

handles.opt.theta = zeros(1, handles.opt.L);
handles.opt.sx = zeros(1, handles.opt.L);
handles.opt.sy = zeros(1, handles.opt.L);
yvecs = handles.yvecs;

%% Register images
for kk = 2: handles.opt.L
    [newsk, Lambdak, yhatk] = LKvar(yvecs{1}, kk,yvecs{kk},0,1,1,zeros(3), 1, handles.opt);
    handles.opt.theta(kk) = newsk(1);
    handles.opt.sx(kk) = handles.opt.res * newsk(2);
    handles.opt.sy(kk) = handles.opt.res * newsk(3);
end

%% 
handles.opt.theta_init = handles.opt.theta;
handles.opt.sx_init = handles.opt.sx;
handles.opt.sy_init = handles.opt.sy;

handles.opt.M = handles.opt.m * handles.opt.res; % set High Resolution image size
handles.opt.N = handles.opt.n * handles.opt.res; % set High Resolution image size

handles.opt.EstimateRegistration = 1;
handles.opt.sigma =0.0;
ADDREGNOISE = 0;

%% Blurring matrix
%% TODO: accept command line parameters
%% TODO: document cases

o_blur = 1 % hardcodes blur option

switch o_blur
    case 1
        handles.opt.h = 1;
    case 2
        size_h = get(handles.edit_size_h_HR,'String');
        size_h = uint8(str2num(size_h));

        if isempty(size_h)
            size_h=0;
        end

        if size_h <= 0
            warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                'modal');
            set(handles.edit_size_h_HR,'String','3');
            handles.opt.h = fspecial('average');
        else
            handles.opt.h = fspecial('average',double(size_h));
        end

    case 3
        size_h = get(handles.edit_size_h_HR,'String');
        size_h = uint8(str2num(size_h));
        sgm = get(handles.edit_varh_HR,'String');
        sgm = str2num(sgm);

        if isempty(size_h)
            size_h=0;
        end

        if isempty(sgm)
            sgm = 0;
        end

        if (size_h <= 0) | (sgm <= 0)
            warndlg({'Invalid blurring matrix parameter values.',' Using default values '},'Invalid Value',...
                'modal');
            set(handles.edit_size_h_HR,'String','3');
            set(handles.edit_varh_HR,'String','0.5');
            handles.opt.h = fspecial('gaussian');
        else
            handles.opt.h = fspecial('gaussian',double(size_h),sgm);
        end

    case 4
        size_h = get(handles.edit_size_h_HR,'String');
        size_h = str2num(size_h);
        
        if isempty(size_h)
            size_h=0;
        end
        
        if size_h <= 0
            warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                'modal');
            set(handles.edit_size_h_HR,'String','5');
            handles.opt.h = fspecial('disk');
        else
            handles.opt.h = fspecial('disk',size_h);
        end

    case 5
       % save('tempSR/h_aux.mat','h');
        load('tempSR/h_aux.mat');
        h =h/(sum(sum(h)));
        handles.opt.h = h;
        
    case 6
        size_h = get(handles.edit_size_h_HR,'String');
        size_h = uint8(str2num(size_h));
        
        if isempty(size_h)
            size_h=0;
        end
        
        xx=reshape(handles.yvecs{1},handles.opt.m,handles.opt.n);
        xx = imresize(xx, handles.opt.res, 'bicubic');
        
        if size_h <= 0
            warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                'modal');
            set(handles.edit_size_h_HR,'String','3');
            
            [xx, handles.opt.h] = deconvblind(xx,ones(size_h,size_h));
           clear xx
        else
            [xx,  handles.opt.h] = deconvblind(xx,ones(size_h,size_h));
            clear xx
        end    
        
        h=handles.opt.h;
        save('h_deconvblind.mat','h');
end
