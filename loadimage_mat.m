%% abre una imagen .mat o de otro tipo y si cancelas devuelve -1

function [image,im_name]= loadimage_mat( filename,pathname)

    f_image = fullfile(pathname, filename);

[PATHSTR,NAME,EXT,VERSN] = fileparts(filename);

load (f_image);
 es=whos ('-file', f_image);
name_im = es(1).name;
 eval(['im= '  name_im ';']);
image= double(im);
im_name = NAME;
end