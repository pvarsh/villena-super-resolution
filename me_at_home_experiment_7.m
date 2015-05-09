opt.filepath = '/home/cusp/pv629/SR/data/me_at_home/7'
out_basepath = '/home/cusp/pv629/SR/output/villena/me_at_home_exp_7' 


%% No blur. Bicubic.
mkdir(out_basepath, '1')
opt.outpath  = fullfile(out_basepath, '1')
opt.sr_method = 1;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 'None';
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))
pause(1); % pause in order to ensure filenames are different

%%------------------- SAR-TV
mkdir(out_basepath, '2')
opt.outpath  = fullfile(out_basepath, '2')
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))

% lambda .5
mkdir(out_basepath, '3')
opt.outpath  = fullfile(out_basepath, '3')
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))

% lambda 1
mkdir(out_basepath, '4')
opt.outpath  = fullfile(out_basepath, '4')
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 1;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))


%%------------------- L4-SAR
mkdir(out_basepath, '5')
opt.outpath  = fullfile(out_basepath, '5')
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))

% lambda .5
mkdir(out_basepath, '6')
opt.outpath  = fullfile(out_basepath, '6')
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))

%%------------------- Robust SR
mkdir(out_basepath, '7')
opt.outpath  = fullfile(out_basepath, '7')
opt.sr_method = 7;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))

%%------------------- Zomet
mkdir(out_basepath, '8')
opt.outpath  = fullfile(out_basepath, '8')
opt.sr_method = 8;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

mkdir(out_basepath, 'ground_truth')
mkdir(out_basepath, 'lr')
copyfile(fullfile(opt.filepath, 'ground_truth'), fullfile(out_basepath, 'ground_truth'))
copyfile(fullfile(opt.filepath, 'face_01_40x40_00.png'), fullfile(out_basepath, 'lr'))
