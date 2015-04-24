%% ----------------------------------------------------------------
%  Parameters:
%       - filepath: directory with low-resolution files
%       - outpath: directory for output
%       - sr_method [1-8]:
%               - 1: Bicubic (using 1st image only)
%               - 2: variational
%               - 3: variational with L4 prior
%               - 4: variational with SAR prior
%               - 5: variational with L4 SAR combination prior
%               - 6: variational TV SAR combination
%               - 7: RobustSR
%               - 8: zomet
%       - blur_method:
%               - 1: No deblurring
%               - 2: Average blur
%               - 3: Gaussian blur
%               - 4: Radial blur
%               - 5: Custom matrix (not tested)
%               - 6: Blind deconvolution (not tested)
%       - varargin:
%               - blur_size (in pixels (TODO: SR or LR pixels?))
%               - blur_variance 
%               - prior distribution weight for mixture of distributions
%                 (e.g. L4 and SAR, TV and SAR etc.) 
%% ----------------------------------------------------------------






% function pipe2(filepath, outpath, sr_method, blur_method, blur_size, blur_var)
% function pipe2(filepath, outpath, sr_method, blur_method, varargin)
function pipe2(pipe_options)
    
    filepath = pipe_options.filepath;
    outpath = pipe_options.outpath;
    handles.opt.outpath = outpath;
    sr_method = pipe_options.sr_method;
    blur_method = pipe_options.blur_method;
    blur_size = pipe_options.blur_size;
    blur_var = pipe_options.blur_var;
    handles.opt.lambda_prior = pipe_options.lambda_prior;

    timestamp = datestr(now);
    
    handles.opt.out_file_prefix = ['/sr_out_' timestamp];

    handles.opt.WriteImages = 1;

    % if nargin >= 5
    %     blur_size = varargin{1};
    % end

    % if nargin >= 6
    %     blur_var = varargin{2};
    % end

    % if nargin == 7
    %     handles.opt.lambda_prior = varargin{3};
    % else
    %     lambda_prior = []
    % end

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
        im_files = dir(strcat(filepath,'/*.png'));
        num_LR = numel(im_files);
        sprintf('>> Loading %d images from argument path...', numel(im_files))
        for i=1:num_LR
            im_f_name = strcat(filepath, '/', im_files(i).name);
            disp(im_files(i).name)

            % grayscale conversion (Algorithms do not use color information)
            im_info = imfinfo(im_f_name);
            im_color_type = im_info.ColorType;
            im = imread(im_f_name);
            if ~strcmp(im_color_type, 'grayscale')
                im = rgb2gray(im);
            end

            % end: grayscale conversion

            if i == 1 % initialize container
                [h, w, chan] = size(im);
                max_dim = max(h,w);
                images = zeros(max_dim, max_dim, num_LR);
                handles.opt.m = h;
                handles.opt.n = w;
            end % end: if

            % square images with letterbox or pillarbox
            [nrow, ncol] = size(im);
            if nrow > ncol
                box_edge = floor((nrow-ncol)/2)+1;
                square_im = zeros(nrow, nrow);
                square_im(:, box_edge:box_edge+ncol-1) = im(:,:);
                im = square_im;
                clear square_im;
            elseif ncol >nrow
                box_edge = floor((ncol-nrow)/2)+1;
                square_im = zeros(ncol, ncol);
                % imshow(square_im)
                square_im(box_edge:box_edge+nrow-1, :) = im(:,:);
                im = square_im;
                clear square_im;
            end
            [handles.opt.m, handles.opt.n] = size(im);
            % end: square images with letterbox or pillarbox

            images(:,:,i) = im;
        end % end: for

    % else % load 4 hardcoded images
    %     disp('>> Loading hardcoded images...')
    %     num_LR = 4;
    %     for i = 1:num_LR
    %         im_f_name = strcat('images/lena1_',num2str(i),'.png');
    %         im = imread(im_f_name);
    %         if not (strcmp())
    %         if i == 1 % initialize container
    %             [h,w,chan] = size(im);
    %             images = zeros(h,w,num_LR);
    %             images(:,:,i) = im;
    %             handles.opt.m = h;
    %             handles.opt.n = w;
    %         end
    %         images(:,:,i) = im;
    %     end % for
    end % if

    y = []; % all lr_images as column vectors

    for k = 1:num_LR
        yk = double(images(:,:,k));
        ys = yk/(max(max(yk))); % normalize

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
    handles.opt.maxit = 40;%opt.maxit = 100;
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
    handles.opt.maxit = 40;
    handles.opt.InitImgExists = 0;
    handles.opt.ShowImages = 0; % was 1
    % handles.opt.WriteImages = 0;
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
    % handles.opt.lambda_prior = 1;
    handles.opt.method = 'variational';

    handles.opt.ApproxSigma = 1;

    handles.opt.DIVIDE_U = 0;
    handles.opt.FixedParameters = 0;
    %% end from LoadParametersSR.m

    % handles.opt.lambda_prior = num2str(0.5); % weight for combination of priors

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

    % blur_method = 2; % hardcodes blur option
    % blur_size = 3;

    % disp('>> Deblurring...');
    switch blur_method
        case 1
            disp('>> No deblurring')
            handles.opt.h = 1;
        case 2
            % blur_size = get(handles.edit_blur_size_HR,'String');
            % blur_size = uint8(str2num(blur_size));
            disp(['>> Deblurring: Average blur with size ' num2str(blur_size)])
            if isempty(blur_size)
                blur_size=0;
            end

            if blur_size <= 0
                warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                    'modal');
                set(handles.edit_blur_size_HR,'String','3');
                handles.opt.h = fspecial('average');
            else
                handles.opt.h = fspecial('average',double(blur_size));
            end

        case 3
            disp(['>> Deblurring: Gaussian blur with size ' num2str(blur_size) ' and variance ' num2str(blur_var)])

            % blur_size = get(handles.edit_blur_size_HR,'String');
            % blur_size = uint8(str2num(blur_size));
            % blur_var = get(handles.edit_varh_HR,'String');
            % blur_var = str2num(blur_var);

            if isempty(blur_size)
                blur_size=0;
            end

            if isempty(blur_var)
                blur_var = 0;
            end

            if (blur_size <= 0) | (blur_var <= 0)
                warndlg({'Invalid blurring matrix parameter values.',' Using default values '},'Invalid Value',...
                    'modal');
                set(handles.edit_blur_size_HR,'String','3');
                set(handles.edit_varh_HR,'String','0.5');
                handles.opt.h = fspecial('gaussian');
            else
                handles.opt.h = fspecial('gaussian',double(blur_size), blur_var);
            end

        case 4
            disp(['>> Deblurring: Radial blur with size ' num2str(blur_size)])
            % blur_size = get(handles.edit_blur_size_HR,'String');
            % blur_size = str2num(blur_size);
            
            if isempty(blur_size)
                blur_size = 0;
            end
            
            if blur_size <= 0
                warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                    'modal');
                set(handles.edit_blur_size_HR,'String','5');
                handles.opt.h = fspecial('disk');
            else
                handles.opt.h = fspecial('disk',blur_size);
            end

        case 5
           % save('tempSR/h_aux.mat','h');
            load('tempSR/h_aux.mat');
            h =h/(sum(sum(h)));
            handles.opt.h = h;
            
        case 6
            blur_size = get(handles.edit_blur_size_HR,'String');
            blur_size = uint8(str2num(blur_size));
            
            if isempty(blur_size)
                blur_size=0;
            end
            
            xx = reshape(handles.yvecs{1},handles.opt.m,handles.opt.n);
            xx = imresize(xx, handles.opt.res, 'bicubic');
            
            if blur_size <= 0
                warndlg({'Invalid blurring matrix size.',' Using default value '},'Invalid Value',...
                    'modal');
                set(handles.edit_blur_size_HR,'String','3');
                
                [xx, handles.opt.h] = deconvblind(xx,ones(blur_size,blur_size));
               clear xx
            else
                [xx,  handles.opt.h] = deconvblind(xx,ones(blur_size,blur_size));
                clear xx
            end    
            
            h=handles.opt.h;
            save('h_deconvblind.mat','h');
    end % end cases

    %% Skipping lines 770-813

    method = sr_method; % hardcodes solvex_var
    % if ~isempty(lambda_prior)
    %     handles.opt.lambda_prior = lambda_prior;
    % end

    % handles.HRimage = axes();
    handles.HRimage = 5;

    switch method
        case 1 %% Bicubic
            disp('>> Superresolving bicubic...')
            [handles.srimage.x,handles.srimage.out] = BicubicSR(handles.y,handles.opt);
            
        case 2 %% solvex_var  
            disp('>> Superresolving TV...')
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_TV_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= solvex_var(handles.y,handles.opt,handles.HRimage,handles);    
                   
        case 3 %% solvex_varL4
            disp('>> Superresolving TV L4...')
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_L4_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= solvex_varL4(handles.y,handles.opt,handles.HRimage,handles); 
           
        case 4 %% solvex_varL4SAR
            disp(['>> Superresolving variational L4 SAR with lambda ' num2str(handles.opt.lambda_prior) '...'])
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_SAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            % handles.opt.lambda_prior = 0; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varL4SAR(handles.y,handles.opt,handles.HRimage,handles);   

        case 5 %% solvex_varL4SAR Combination (same as 4?)
            disp(['>> Superresolving variational L4 SAR combination with lambda ' num2str(handles.opt.lambda_prior) '...'])         
            FILE_log = fopen(sprintf('tempSR/LOG_VAR_combL4SAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', handles.opt.maxit, handles.opt.PCG_maxit, handles.opt.LK_maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            % handles.opt.lambda_prior = str2num(get(handles.edit_lambda,'string'));

            % set(handles.edit_lambda,'String',num2str(handles.opt.lambda_prior));
            % handles.opt.lambda_prior = 0.5; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varL4SAR(handles.y,handles.opt,handles.HRimage,handles);   
                      
        case 6 %% solvex_varTVSAR Combination
            disp(['>> Superresolving variational TV SAR combination with lambda ' num2str(handles.opt.lambda_prior) '...'])
            FILE_log = fopen(strcat(handles.opt.outpath, ...
                                    handles.opt.out_file_prefix, ...
                                    sprintf('LOG_VAR_combTVSAR_maxit%d_PCGmaxit%d_LKmaxit%d_sigma%g_RegERR%d_exp1.txt', ...
                                            handles.opt.maxit, ...
                                            handles.opt.PCG_maxit, ...
                                            handles.opt.LK_maxit, ...
                                            handles.opt.sigma, ...
                                            ADDREGNOISE ...
                                            ) ...
                                    ), ...
                            'w' ...
                            );
            handles.opt.LogFile = FILE_log;
            % handles.opt.lambda_prior = str2num(get(handles.edit_lambda,'string'));
            % set(handles.edit_lambda,'String',num2str(handles.opt.lambda_prior));
            % handles.opt.lambda_prior = 0.5; % SAR prior is seleted; with value 1 the norm l1 prior is selected
            [handles.srimage.x,handles.srimage.out]= solvex_varTVSAR(handles.y,handles.opt,handles.HRimage,handles);   
                            
        case 7 %% RobustSR
            disp('>> Superresolving ROBUST SR...')
            FILE_log = fopen(sprintf('tempSR/LOG_RSR_maxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= RobustSR(handles.y,handles.opt,handles.HRimage,handles);
           
        case 8 %% Zomet
            disp('>> Superresolving Zomet...')
            FILE_log = fopen(sprintf('tempSR/LOG_Zomet_maxit%d_sigma%g_RegERR%d.txt', handles.opt.maxit, handles.opt.sigma,ADDREGNOISE),'w');
            handles.opt.LogFile = FILE_log;
            [handles.srimage.x,handles.srimage.out]= Zomet(handles.y,handles.opt,handles.HRimage,handles);          
    end % end switch

    % Create list of SR and blur option codes
    % Note: vertcat function is implicitly used. Requires strings to be of same length
    sr_methods = cellstr(['bicubic        '; ...
                          'sovlex_var     '; ...
                          'solvex_varL4   '; ...
                          'solvex_varL4SAR'; ...
                          'solvex_varL4SAR'; ...
                          'solvex_varTVSAR'; ...
                          'RobustSR       '; ...
                          'Zomet          '; ...
                          ]);
    blur_methods = cellstr(['no blur      '; ...
                            'average blur '; ...
                            'gaussian blur'; ...
                            'radial blur  '; ...
                            ]);

    disp('>> Saving SR image and writing log');
    disp([outpath handles.opt.out_file_prefix '.log']);
    logFileId = fopen([outpath handles.opt.out_file_prefix '.log'], 'w');
    fprintf(logFileId, 'Villena et al. Super Resolution Software. File: sr_cli.m\n');
    fprintf(logFileId, strcat('-', timestamp, '\n'));
    fprintf(logFileId, strcat('-', 'filepath:\t', filepath, '\n'));
    fprintf(logFileId, strcat('-', 'outpath:\t', outpath, '\n'));
    fprintf(logFileId, strcat('-', 'SR filename:\t', 'sr_out_', timestamp, '.png', '\n'));
    fprintf(logFileId, strcat('-', '# images:\t', num2str(handles.opt.L), '\n'));
    fprintf(logFileId, strcat('-', 'srmethod:\t', num2str(sr_method), ' ', sr_methods{sr_method}, '\n'));
    fprintf(logFileId, strcat('-', 'blur:\t\t', num2str(blur_method), ' ', blur_methods{blur_method}, '\n'));
    if ~exist('blur_size')
        blur_size = 'None';
    end
    if ~exist('blur_var')
        blur_var = 'None';
    end
    % disp(class(handles.opt))
    % disp(handles.opt.lambda_prior)
    if ~isfield(handles.opt, 'lambda_prior')
        handles.opt.lambda_prior = 'None';
    end
    fprintf(logFileId, strcat('-', 'blur_size:\t', num2str(blur_size), '\n'));
    fprintf(logFileId, strcat('-', 'blur_variance:\t', num2str(blur_var), '\n'));
    fprintf(logFileId, strcat('-', 'lambda_prior:\t', num2str(handles.opt.lambda_prior), '\n'));

    % fprintf(logFileId, strcat('varargin:', varargin))

    imFileName = [outpath '/sr_out_' timestamp '.png'];

    x = handles.srimage.x;
    x = reshape(x, handles.opt.M, handles.opt.N);
    x = uint8(x);
    imwrite(x, imFileName);

    % close(handles.HRimage)
    %% Show SR image
    % figure('name', 'Myfig', 'Position', [100,100,100,100], 'units', 'normalized', 'outerposition', [0 0 1 1]), imshow(x)  % subplots_adjust(left=None, bottom=None, right=None, top=None,
    %                   wspace=None, hspace=None);
    % figure, imshow(x, 'InitialMagnification', 400);
end % end pipe2
