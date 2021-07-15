%% read images and downsample them
%there used to be a subfolder within the data folder called exposure_stack, 
%where exist 16 RAW files, 16 JPEG files, and 16 TIFF files output from doing 
%dcraw in terminal, these are removed for the sake that the whole folder is too bulky
%That is why the path is shown as below
expo1 = double(imread('../data/Part4/photo1.tiff'));
expo2 = double(imread('../data/Part4/photo2.tiff'));
expo3 = double(imread('../data/Part4/photo3.tiff'));
expo4 = double(imread('../data/Part4/photo4.tiff'));
expo5 = double(imread('../data/Part4/photo5.tiff'));
expo6 = double(imread('../data/Part4/photo6.tiff'));
expo7 = double(imread('../data/Part4/photo7.tiff'));
expo8 = double(imread('../data/Part4/photo8.tiff'));
expo9 = double(imread('../data/Part4/photo9.tiff'));
expo10 = double(imread('../data/Part4/photo10.tiff'));
expo11 = double(imread('../data/Part4/photo11.tiff'));
expo12 = double(imread('../data/Part4/photo12.tiff'));
%%
exposure1 = [expo1(:,:,1),expo1(:,:,2),expo1(:,:,3)];
exposure2 = [expo2(:,:,1),expo2(:,:,2),expo2(:,:,3)];
exposure3 = [expo3(:,:,1),expo3(:,:,2),expo3(:,:,3)];
exposure4 = [expo4(:,:,1),expo4(:,:,2),expo4(:,:,3)];
exposure5 = [expo5(:,:,1),expo5(:,:,2),expo5(:,:,3)];
exposure6 = [expo6(:,:,1),expo6(:,:,2),expo6(:,:,3)];
exposure7 = [expo7(:,:,1),expo7(:,:,2),expo7(:,:,3)];
exposure8 = [expo8(:,:,1),expo8(:,:,2),expo8(:,:,3)];
exposure9 = [expo9(:,:,1),expo9(:,:,2),expo9(:,:,3)];
exposure10 = [expo10(:,:,1),expo10(:,:,2),expo10(:,:,3)];
exposure11 = [expo11(:,:,1),expo11(:,:,2),expo11(:,:,3)];
exposure12 = [expo12(:,:,1),expo12(:,:,2),expo12(:,:,3)];


%% RAW linear uniform
zmin = 655.36;
zmax = 64880.64;
B = 1/2048 * [2^0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11];
%%
[urw1 num1] = uniweight(exposure1, zmin, zmax,B(1));
[urw2 num2] = uniweight(exposure2, zmin, zmax,B(2));
[urw3 num3] = uniweight(exposure3, zmin, zmax,B(3));
[urw4 num4] = uniweight(exposure4, zmin, zmax,B(4));
%%
[urw5 num5] = uniweight(exposure5, zmin, zmax,B(5));
[urw6 num6] = uniweight(exposure6, zmin, zmax,B(6));
[urw7 num7] = uniweight(exposure7, zmin, zmax,B(7));
[urw8 num8] = uniweight(exposure8, zmin, zmax,B(8));
%%
[urw9 num9] = uniweight(exposure9, zmin, zmax,B(9));
[urw10 num10] = uniweight(exposure10, zmin, zmax,B(10));
[urw11 num11] = uniweight(exposure11, zmin, zmax,B(11));
[urw12 num12] = uniweight(exposure12, zmin, zmax,B(12));

%%
numerator = num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12;
denominator = urw1 + urw2 + urw3 + urw4 + urw5 + urw6 + urw7 + urw8 + urw9 + urw10 + urw11 + urw12;

%% RAW Log photon
Myimg = exp(numerator ./ denominator);
save('My_img.mat','Myimg');

%% RAW Log part5
part5 = exp(numerator ./ denominator);
save('part5.mat','part5');

%% Evaluation
% rgb(:,:,1) = I_raw_l_uni(:,1:6016)/255;
% rgb(:,:,2) = I_raw_l_uni(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_l_uni(:,12033:end)/255;
% 
% rgb(:,:,1) = I_raw_l_t(:,1:6016)/255;
% rgb(:,:,2) = I_raw_l_t(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_l_t(:,12033:end)/255;

% rgb(:,:,1) = I_raw_l_g(:,1:6016)/255;
% rgb(:,:,2) = I_raw_l_g(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_l_g(:,12033:end)/255;
%  
% rgb(:,:,1) = I_raw_l_p(:,1:6016)/255;
% rgb(:,:,2) = I_raw_l_p(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_l_p(:,12033:end)/255;

% rgb(:,:,1) = I_raw_lo_u(:,1:6016)/255;
% rgb(:,:,2) = I_raw_lo_u(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_lo_u(:,12033:end)/255;
% 
% rgb(:,:,1) = I_raw_lo_t(:,1:6016)/255;
% rgb(:,:,2) = I_raw_lo_t(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_lo_t(:,12033:end)/255;

% rgb(:,:,1) = I_raw_lo_g(:,1:6016)/255;
% rgb(:,:,2) = I_raw_lo_g(:,6017:12032)/255;
% rgb(:,:,3) = I_raw_lo_g(:,12033:end)/255;
% 
rgb(:,:,1) = part5(:,1:6016)/255;
rgb(:,:,2) = part5(:,6017:12032)/255;
rgb(:,:,3) = part5(:,12033:end)/255;



%% save hdr file

index = find(isnan(rgb) == 1);
rgb(index) = 1;
hdrwrite(rgb/10000,'part5.hdr');
% hdrwrite(rgb/100000,'myimg1.hdr');

 
%% Tone Mapping RGB as a whole
maxi = max(max(max(rgb)));
kkk = rgb/maxi;
K = 0.0005;
B = 0.95;
N = size(kkk,1) * size(kkk,2);
inter0 = kkk(:,:,1) + 1e-200; %r
inter2 = kkk(:,:,2) + 1e-200; %g
inter3 = kkk(:,:,3) + 1e-200; %b
inter1 = [inter0,inter2,inter3];

ImHDR = exp(1/N * sum(sum(log(inter1))));
I_ijHDR = K/ImHDR * inter1;
I_white = B * max(max(I_ijHDR));
I_ijTM = I_ijHDR .* (1+I_ijHDR/I_white/I_white)./(1+I_ijHDR);

FINALIMG(:,:,1) = I_ijTM(:,1:6016);
FINALIMG(:,:,2) = I_ijTM(:,6017:12032);
FINALIMG(:,:,3) = I_ijTM(:,12033:end);
imshow(FINALIMG);
imwrite(FINALIMG,'rgbmyimg000595.png');

%% Tone mapping with Y
maxi = max(max(max(rgb)));
kkk = rgb/maxi;
K = 0.5;
B = 950;
N = size(kkk,1) * size(kkk,2);

zzz = rgb2xyz(kkk,'ColorSpace','linear-rgb');
x = zzz(:,:,1)./(zzz(:,:,1)+zzz(:,:,2)+zzz(:,:,3));
y = zzz(:,:,2)./(zzz(:,:,1)+zzz(:,:,2)+zzz(:,:,3));
Y = zzz(:,:,2);


ImHDR = exp(1/N * sum(sum(log(Y))));
I_ijHDR = K/ImHDR * Y;
I_white = B * max(max(I_ijHDR));
I_ijTM = I_ijHDR .* (1+I_ijHDR/I_white/I_white)./(1+I_ijHDR);


[Xo,Yo,Zo] = xyY_to_XYZ( x, y, I_ijTM );
IIII(:,:,1) = Xo;
IIII(:,:,2) = Yo;
IIII(:,:,3) = Zo;
FINALIMG = xyz2rgb(IIII,'ColorSpace','linear-rgb');

imshow(FINALIMG);
imwrite(FINALIMG,'myimgY05950.png');


