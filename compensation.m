%contrast enhancement
function [result_image blurim]=compensation(input_image,padlength,op_cen_out,op_cen_in,...
    op_sur_out,op_sur_in,op_rmt_out,...
    op_rmt_in);

%Expand pic to void edge problems:
input_image=expandpic(input_image,padlength);                          

%Extract RGB value's from input image
R=double(input_image(:,:,1)); 
G=double(input_image(:,:,2)); 
B=double(input_image(:,:,3));

%Transformation from RGB to L which is the luminance (Y in XYZ)
L=(0.177*R+0.813*G+0.011*B);
% L = R+G+B;
L(find(L==0))=1e-6;       %this is to void division in zero later

LRGB=255;
Rcone=R./LRGB;
Gcone=G./LRGB;
Bcone=B./LRGB;
Ucone=L./LRGB;         %Chromatic Channel


%Definitions OPPONENT RECEPTIVE FIELDS 
%--------------------------------------------------------------------
center=rfields(op_cen_out,op_cen_in);    %center mask area
surround=rfields(op_sur_out,op_sur_in);  %the surround area mask 5,1
remote=rfields(op_rmt_out,op_rmt_in);    %remote area mask 30,3

blue_center=blue_fields(15*op_cen_out,8*op_cen_out);
center_blur=rfields(5*op_cen_out,op_cen_in);%9    %center mask area

%--------------------------------------------------------------------


%CREATE OPONNENT SIGNAL
%-----------------------------------------------
Rcenter=filter2(center,Rcone);
Gcenter=filter2(center,Gcone);

Bcenter=filter2(center,Bcone);

 
Mcenter = (Rcenter+Bcenter)/2;
Mbcenter = filter2(blue_center,Mcenter);
Bbcenter = filter2(blue_center,Bcenter);
Gbcenter = filter2(blue_center,Gcenter);
Rbcenter = filter2(blue_center,Rcenter);

Ubcenter = filter2(blue_center,Ucone);

Rs=filter2(surround,Rcone);
Gs=filter2(surround,Gcone);   % green surround map

Us=filter2(surround,Ucone);

%SIGNAL AFTER ADAPTATION:
A = 1;
C_w = 1;
C_remote = 1;
Center_Surround_Ratio = 1;
Center_Surround_Total = 2;
BB = Center_Surround_Total/(1+Center_Surround_Ratio);
AA = Center_Surround_Total - B;

AA = 1;
BB = 1;

RGon = Rcenter-Gs;
GRon = Gcenter-Rs;
BYon = Bbcenter -(Rbcenter+Gbcenter)/2;


Rnew = (Rcenter+Gcenter)/2;
Gnew = Rnew;
Bnew = Rnew;
RGon= Rcenter-Gs;
GRon= Gcenter-Rs;
BYon= Bbcenter-(Rbcenter+Gbcenter)/2;


%Jacobi Diffusion
for(i = 1:100)
    Rs=filter2(surround,Rnew);
    Gs=filter2(surround,Gnew);   % green surround map
    
    Rbcenter = filter2(blue_center,Rnew);
    Gbcenter = filter2(blue_center,Gnew);
    Bbcenter=filter2(blue_center,Bcenter);
    Ybcenter = (Gbcenter+Rbcenter)/2;
    
    dt = 0.05;
    Rnew = Rnew-dt*(Rnew-(RGon+Gs));
    Gnew = Gnew-dt*(Gnew-(GRon+Rs));
    Bnew = BYon+(Rnew+Gnew)/2;
    saver(i)= Rnew(60,60);
    err2(i)= (Rnew(60,60)-(RGon(60,60)+Gs(60,60)));
end

im0(:,:,1)=Rnew.*LRGB;
im0(:,:,2)=Gnew.*LRGB;
im0(:,:,3)=Bnew.*LRGB;

Ln=(0.177*im0(:,:,1)+0.813*im0(:,:,2)+0.011*im0(:,:,3));
for (i=1:1:3) im0(:,:,i)=im0(:,:,i).*L./Ln; end

im0=de_expand(im0,padlength);
% figure(2)
blurim(:,:,1)=Rcenter.*LRGB;
blurim(:,:,2)=Gcenter.*LRGB;
blurim(:,:,3)=Bbcenter.*LRGB;
blurim=de_expand(blurim,padlength);
blurim = uint8(blurim);
result_image=uint8(im0);

end

function expandpic_result=expandpic(im,d)
[xmax,ymax]=size(im(:,:,1));
for(i=-d+1:1:(xmax+d))
    for(j=-d+1:1:(ymax+d))
        if i<1 in=2-i;
        elseif i>xmax in=2*xmax-i;
        else in=i;
        end
        if j<1 jn=2-j;
        elseif j>ymax jn=2*ymax-j;
        else jn=j;
        end
        imn(i+d,j+d,:)=im(in,jn,:);
    end
end
expandpic_result=imn;
end

function de_expand_result=de_expand(im,d)
[m,n]=size(im(:,:,1));
de_expand_result(1:1:m-2*d,1:1:n-2*d,:)=im(1+d:1:m-d,1+d:1:n-d,:);
end

%this function returns a sizeXsize matrix that contains the receptive field mask , using a gaussian
%inner size is the size of matrix that will be zero insize the sizeXsize matrix
%for example by choosing size=3, isize=1, Kfield=infinity function will return 1/8[1 1 1;1 0 1;1 1 1];
%the function is also normalizing the matrix

function rfields_result=rfields(size,isize,Kfield);
field=zeros(size,size);
%Kfield exponent decaying constant
if isize~=0
    for(i=1:1:size)
     for(j=1:1:size)
           field(i,j)=-exp(-( (i-(size+1)/2)^2+(j-(size+1)/2)^2)/((isize+0.0001)^2))+exp(-( (i-(size+1)/2)^2+(j-(size+1)/2)^2)/(0.1*(size+0.001)^2));    %mask is givven the value of decaying
      
      end                                                    %exponent according to the distance from  
    end
else
     for(i=1:1:size)
     for(j=1:1:size)
           field(i,j)=exp(-( (i-(size+1)/2)^2+(j-(size+1)/2)^2)/((size)^2));    %mask is givven the value of decaying
      
      end                                                    %exponent according to the distance from  center
    end
end

field(round((size-isize+2)/2):round((size+isize)/2),round((size-isize+2)/2):round((size+isize)/2))=zeros(isize,isize);         %zero padding in innersize
field=field/(sum(sum(field)));  %normalization to 1
rfields_result=field;
end

function blue_fields_result=blue_fields(size,Kfield);
field=zeros(size,size);
index = (1:size) - (1+size)/2;
[i j] = meshgrid(index,index);

field =exp(-(i.^2+j.^2)./Kfield.^2);    %mask is givven the value of decaying

field=field/(sum(sum(field)));  %normalization to 1
blue_fields_result=field;
end

