function pipe2(filepath, outpath, sr_method)

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
    %% TODO: grayscale conversion

    %% Load files from directory
    if (filepath ~= ' ') % load from directory
        disp('>> Loading images from argument path...')
        im_files = dir(strcat(filepath,'/*.png'));
        num_LR = numel(im_files);
        for i=1:num_LR
            im_f_name = strcat(filepath, '/', im_files(i).name);
            im = imread(im_f_name);
            if i == 1 % initialize container
                [h,w,chan] = size(im);
                images = zeros(h,w,num_LR);
                images(:,:,i) = im;
                handles.opt.m = h;
                handles.opt.n = w;
            end % if
        end % for
            images(:,:,i) = im;
    else % load 4 hardcoded images
        disp('>> Loading hardcoded images...')
        num_LR = 4;
        for i = 1:num_LR
            im_f_name = strcat('images/lena1_',num2str(i),'.png');
            im = imread(im_f_name);
            if i == 1 % initialize container
                [h,w,chan] = size(im);
                images = zeros(h,w,num_LR);
                images(:,:,i) = im;
                handles.opt.m = h;
                handles.opt.n = w;
            end
            images(:,:,i) = im;
        end % for
    end % if

    y = []; % all lr_images as column vectors

    for k = 1:num_LR
        yk = double(images(:,:,k));
        % ys = yk/(max(max(yk))); % normalize

        yvecs{k} = yk(:);
        y = [y; yk(:)];
        handles.y = y;
        handles.opt.y = y;
    end

    handles.yvecs = yvecs;

    handles.opt.L = num_LR;
    handles.opt.res = 2;

    %% Pre-LK registration options

    % Load Base Options
    %% Specific to Bayesian methods
    opt.noise_estimate = 'COMMON';
    handles.opt.noise_estimate = 'SEPARATE';
    handles.opt.LK_maxit = 35;
    handles.opt.LK_thr = 1e-4;
    handles.opt.PCG_minit = 20;
    handles.opt.PCG_maxit = 100;
    handles.opt.PCG_thr = 1e-10;
    handles.opt.maxit = 100;%opt.maxit = 100;
    % handles.opt.thr = 1e-5;

    handles.opt.InitWithAvgImg = 1;
    handles.opt.InitImgExists = 0;

    % %% Output options
    % handles.opt.WriteImages = 0;
    % handles.opt.ShowImages = 1;
    % handles.opt.KeepHistory = 0;
    % handles.opt.verbose = 1; 
    handles.opt.debug = 0;
    handles.opt.interv_ShowImages = 0; %was 1 % Show images each opt.interv_ShowImages iterations
    handles.opt.FixedBeta = 0;        
    %----------------

    %% From LoadParametersSR.m
    % handles.opt.P = MFSR;
    % handles.output = hObject;
    % guidata(hObject, handles);
    % handles.opt.sx = handles.opt.sx_init;
    % handles.opt.sy = handles.opt.sy_init;
    % handles.opt.theta = handles.opt.theta_init;
    handles.opt.P = double(handles.opt.res);
    handles.opt.alpha = 0.6;
    handles.opt.beta = 1/100;
    handles.opt.lambda = 0.005;
    handles.opt.maxit = 300;
    handles.opt.InitImgExists = 0;
    handles.opt.ShowImages = 0; % was 1
    handles.opt.WriteImages = 0;
    handles.opt.KeepHistory = 0;
    handles.opt.thr = 1e-5;
    handles.opt.PCG_maxit = 100;
    handles.opt.LK_maxit = 100;
    handles.opt.debug = 0;
    handles.opt.LK_thr =1e-4;
    handles.opt.noise_estimate = 'SEPARATE';
    handles.opt.PCG_thr =1e-10;
    handles.opt.PCG_minit = 10;
    handles.opt.verbose = 1;
    handles.opt.FixedBeta = 0 ; 
    handles.opt.lambda_prior = 1;
    handles.opt.method = 'variational';

    handles.opt.ApproxSigma = 1;

    handles.opt.DIVIDE_U = 0;
    handles.opt.FixedParameters = 0;
    %% end from LoadParametersSR.m

    handles.opt.lambda_prior = num2str(0.5); % weight for combination of priors

    % set(handles.text_Runing,'Visible','Off');
    % set(handles.text_Runing,'Enable','Off');
    % set(handles.opt.ref_cancel,'Value',0);
    handles.opt.ref_cancel = 0;
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

    disp('>> Registering...')
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

    o_blur = 1; % hardcodes blur option

    disp('>> Deblurring...');
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
    end % end cases

    %% Skipping lines 770-813

    method = sr_method; % hardcodes solvex_var

    handles.HRimage = axes();

    disp('>> Superresolving...')
    switch method
        case 1 %% Bicubic
            [handles.srimage.x,handles.srimage.out] = BicubicSR(handles.y,handles.opt);
            
        case 2 %% solvex_var    
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_TV_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= solvex_var(handles.y,handles.opt,handles.HRimage,handles);    
                   
        case 3 %% solvex_varL4
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_L4_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= solvex_varL4(handles.y,handles.opt,handles.HRimage,handles); 
           
        case 4 %% solvex_varL4SAR
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_SAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            handles.opt.lambda_prior = 0; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varL4SAR(handles.y,handles.opt,handles.HRimage,handles);   

        case 5 %% solvex_varL4SAR Combination            
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_combL4SAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            % handles.opt.lambda_prior = str2num(get(handles.edit_lambda,'string'));

            % set(handles.edit_lambda,'String',num2str(handles.opt.lambda_prior));
            %handles.opt.lambda_prior = 0.5; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varL4SAR(handles.y,handles.opt,handles.HRimage,handles);   
                      
        case 6 %% solvex_varTVSAR Combination
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_combTVSAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            % handles.opt.lambda_prior = str2num(get(handles.edit_lambda,'string'));
            % set(handles.edit_lambda,'String',num2str(handles.opt.lambda_prior));
            %handles.opt.lambda_prior = 0.5; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varTVSAR(handles.y,handles.opt,handles.HRimage,handles);   
                            
        case 7 %% RobustSR
            FILE_log = fopen(sprintf('tempSR/LOG_RSR_maxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= RobustSR(handles.y,handles.opt,handles.HRimage,handles);
           
        case 8 %% Zomet           
            FILE_log = fopen(sprintf('tempSR/LOG_Zomet_maxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= Zomet(handles.y,handles.opt,handles.HRimage,handles);          
    end % end switch

    disp('>> Saving SR image')
    x = handles.srimage.x;
    x = reshape(x, handles.opt.M, handles.opt.N);
    x = uint8(x);
    imwrite(x, strcat(outpath, '/out_', datestr(now), '.png'))

    %% Show SR image
    figure(), imshow(x);

end % end pipe2
