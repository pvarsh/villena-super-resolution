

inpath = '/Users/petervarshavsky/Dropbox/NYU/superresolution/data/CCTV-robbery-youtube/extracted_sequences/license_plate'
outpath = '/Users/petervarshavsky/Dropbox/NYU/superresolution/output/villena/CCTV_lp_0321'

% blurs = {1,2,3,4}
% sr_methods = {1,2,3,4}

blur = 2
blur_size = 3
sr_method = 2

sr_cli(inpath, outpath, sr_method, blur, blur_size)