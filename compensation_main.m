function compensation_main()
%load image

name = 'oi_image.tif'

close all;
out_file=[name(1:end-4) '_fix' '.bmp'];  %output filename
blur_file=[name(1:end-4) '_blur' '.bmp'];  %output filename

filetype=name(end-2:end);                                         %file type can be almost any type (such as bmp)
input_image=imread(name,filetype);                %reads the image to a 3D matrix consist of RGB value

figure(1)
image(input_image);
title('Input Image');
                                                    %for any pixel
                                                    
%parameters:
padlength=60;
op_cen_out=1;
op_cen_in=0;
op_sur_out=3;
op_sur_in=1;
op_rmt_out=35;
op_rmt_in=3;


[im0 blurim]=compensation(input_image,padlength,...
    op_cen_out,op_cen_in,op_sur_out,op_sur_in,op_rmt_out,op_rmt_in);

figure (2)
title('Output Image');
image(uint8(im0));


imwrite(uint8(im0),out_file);
imwrite(uint8(blurim),blur_file);


