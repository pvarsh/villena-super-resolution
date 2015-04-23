opt.filepath = '/home/cusp/pv629/SR/data/CCTV-robbery-youtube/extracted_sequences/license_plate_80x80_long'
opt.outpath = '/home/cusp/pv629/SR/output/villena/lp_80x80_long'

%% No blur. Bicubic.
opt.sr_method = 1;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 'None';
sr_cli(opt);

%%------------------- SAR-TV
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = 0;
sr_cli(opt);

% lambda .5
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = .5;
sr_cli(opt);

% lambda .5

% lambda 1
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 6;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = 1;
sr_cli(opt);



%%------------------- L4-Sar
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = 0;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = 0;
sr_cli(opt);

% lambda .5
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = .5;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = .5;
sr_cli(opt);

% lambda .5

% lambda 1
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 1;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 3;
opt.blur_var = 3;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 1;
opt.lambda_prior = 1;
sr_cli(opt);

opt.sr_method = 5;
opt.blur_method = 2;
opt.blur_size = 5;
opt.blur_var = 3;
opt.lambda_prior = 1;
sr_cli(opt);


