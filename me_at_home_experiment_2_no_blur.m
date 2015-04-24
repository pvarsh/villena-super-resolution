opt.filepath = '/home/cusp/pv629/SR/data/me_at_home/2'
opt.outpath = '/home/cusp/pv629/SR/output/villena/me_at_home_exp_2'

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

% lambda .5
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);


% lambda 1
opt.sr_method = 6;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 1;
sr_cli(opt);



%%------------------- L4-SAR
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);


% lambda .5
opt.sr_method = 5;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = .5;
sr_cli(opt);

%%------------------- Robust SR
opt.sr_method = 7;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

%%------------------- Zomet
opt.sr_method = 8;
opt.blur_method = 1;
opt.blur_size = 'None';
opt.blur_var = 'None';
opt.lambda_prior = 0;
sr_cli(opt);

