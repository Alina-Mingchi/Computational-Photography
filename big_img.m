%% read images and downsample them
%there used to be a subfolder within the data folder called exposure_stack, 
%where exist 16 RAW files, 16 JPEG files, and 16 TIFF files output from doing 
%dcraw in terminal, these are removed for the sake that the whole folder is too bulky
%That is why the path is shown as below
expo1 = double(imread('../data/exposure_tiff/exposure1.tiff'));
expo2 = double(imread('../data/exposure_tiff/exposure2.tiff'));
expo3 = double(imread('../data/exposure_tiff/exposure3.tiff'));
expo4 = double(imread('../data/exposure_tiff/exposure4.tiff'));
expo5 = double(imread('../data/exposure_tiff/exposure5.tiff'));
expo6 = double(imread('../data/exposure_tiff/exposure6.tiff'));
expo7 = double(imread('../data/exposure_tiff/exposure7.tiff'));
expo8 = double(imread('../data/exposure_tiff/exposure8.tiff'));
expo9 = double(imread('../data/exposure_tiff/exposure9.tiff'));
expo10 = double(imread('../data/exposure_tiff/exposure10.tiff'));
expo11 = double(imread('../data/exposure_tiff/exposure11.tiff'));
expo12 = double(imread('../data/exposure_tiff/exposure12.tiff'));
expo13 = double(imread('../data/exposure_tiff/exposure13.tiff'));
expo14 = double(imread('../data/exposure_tiff/exposure14.tiff'));
expo15 = double(imread('../data/exposure_tiff/exposure15.tiff'));
expo16 = double(imread('../data/exposure_tiff/exposure16.tiff'));

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
exposure13 = [expo13(:,:,1),expo13(:,:,2),expo13(:,:,3)];
exposure14 = [expo14(:,:,1),expo14(:,:,2),expo14(:,:,3)];
exposure15 = [expo15(:,:,1),expo15(:,:,2),expo15(:,:,3)];
exposure16 = [expo16(:,:,1),expo16(:,:,2),expo16(:,:,3)];


%% RAW linear uniform
zmin = 655.36;
zmax = 64880.64;
B = 1/2048 * [2^0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14,2^15];
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
[urw13 num13] = uniweight(exposure13, zmin, zmax,B(13));
[urw14 num14] = uniweight(exposure14, zmin, zmax,B(14));
[urw15 num15] = uniweight(exposure15, zmin, zmax,B(15));
[urw16 num16] = uniweight(exposure16, zmin, zmax,B(16));
%%
numerator = num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14+num15+num16;
denominator = urw1 + urw2 + urw3 + urw4 + urw5 + urw6 + urw7 + urw8 + urw9 + urw10 + urw11 + urw12 + urw13 + urw14 + urw15 + urw16;
%% RAW Linear Uniform
I_raw_l_uni = numerator ./ denominator;
save('raw_li_uni.mat','I_raw_l_uni');
%% RAW Linear Tent
I_raw_l_t = numerator ./ denominator;
save('raw_li_t.mat','I_raw_l_t');
%% RAW Linear gaussian
I_raw_l_g = numerator ./ denominator;
save('raw_li_g.mat','I_raw_l_g');
%% RAW Linear photon
I_raw_l_p = numerator ./ denominator;
save('raw_li_p.mat','I_raw_l_p');
%% RAW Log uni
I_raw_lo_u = exp(numerator ./ denominator);
save('raw_lo_u.mat','I_raw_lo_u');
%% RAW Log tent
I_raw_lo_t = exp(numerator ./ denominator);
save('raw_lo_t.mat','I_raw_lo_t');
%% RAW Log gaussian
I_raw_lo_g = exp(numerator ./ denominator);
save('raw_lo_g.mat','I_raw_lo_g');
%% RAW Log photon
I_raw_lo_p = exp(numerator ./ denominator);
save('raw_lo_p.mat','I_raw_lo_p');

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
rgb(:,:,1) = I_raw_lo_p(:,1:6016)/255;
rgb(:,:,2) = I_raw_lo_p(:,6017:12032)/255;
rgb(:,:,3) = I_raw_lo_p(:,12033:end)/255;

img_XYZ = rgb2xyz(rgb,'ColorSpace','linear-rgb');

%read color checker and display
[ r, g, b] = read_colorchecker_gm();
checker(:,:,1) = r;
checker(:,:,2) = g;
checker(:,:,3) = b;
% figure
% imshow(checker);
checker_XYZ = rgb2xyz(checker,'ColorSpace','linear-rgb');
% get the square
square = img_XYZ(580:1564,3260:3917,:);
figure
imshow(square);
impixelinfo()
%
ps_y = [mean(mean(square(824:956,512:640,2))),mean(mean(square(684:796,506:636,2))),...
    mean(mean(square(506:636,502:628,2))),mean(mean(square(352:476,498:622,2))),...
    mean(mean(square(194:320,494:616,2))),mean(mean(square(38:162,490:608,2)))];
x = 1:6;
y = log(ps_y);
bb = polyfit(x, y, 1);
%bb(1) slope, bb(2) intercept
yCalc1 = x*bb(1)+bb(2);
scatter(x,y)
hold on
plot(x,yCalc1)
xlabel('number')
ylabel('log value of illuminance')
title('plot of log regression of RAW photon weights')
grid on

error = sum((y - yCalc1).^2);

%% save hdr file

index = find(isnan(rgb) == 1);
rgb(index) = 1;
    
hdrwrite(rgb/100,'before.hdr');

%% color correction
square = rgb(580:1564,3260:3917,:);
figure;
imshow(square);
ch1 = mean(mean(square(828:960,44:168,:)));
ch2 = mean(mean(square(828:958,200:328,:)));
ch3 = mean(mean(square(834:958,362:484,:)));
ch4 = mean(mean(square(824:956,512:640,:)));
ch5 = mean(mean(square(666:796,36:162,:)));
ch6 = mean(mean(square(672:794,198:322,:)));
ch7 = mean(mean(square(666:794,356:478,:)));
ch8 = mean(mean(square(684:796,506:636,:)));
ch9 = mean(mean(square(506:636,34:156,:)));
ch10 = mean(mean(square(506:634,190:316,:)));
ch11 = mean(mean(square(508:634,348:474,:)));
ch12 = mean(mean(square(506:636,502:628,:)));
ch13 = mean(mean(square(350:474,28:154,:)));
ch14 = mean(mean(square(350:474,188:310,:)));
ch15 = mean(mean(square(352:474,340:468,:)));
ch16 = mean(mean(square(352:476,498:622,:)));
ch17 = mean(mean(square(186:316,22:148,:)));
ch18 = mean(mean(square(188:314,180:304,:)));
ch19 = mean(mean(square(192:316,336:460,:)));
ch20 = mean(mean(square(194:320,494:616,:)));
ch21 = mean(mean(square(28:158,16:144,:)));
ch22 = mean(mean(square(30:158,172:300,:)));
ch23 = mean(mean(square(36:160,330:454,:)));
ch24 = mean(mean(square(38:162,490:608,:)));
%%
average_square = [ch1(:),ch2(:),ch3(:),ch4(:),ch5(:),ch6(:),ch7(:),ch8(:),ch9(:),ch10(:),...
    ch11(:),ch12(:),ch13(:),ch14(:),ch15(:),ch16(:),ch17(:),ch18(:),ch19(:),ch20(:),...
    ch21(:),ch22(:),ch23(:),ch24(:)];
real_r = r(:)';
real_g = g(:)';
real_b = b(:)';
%%
term = ones(1,24);
real_checker = [r(:)';g(:)';b(:)';term];
a_s = average_square';
r_c = real_checker';
T = (r_c \ a_s)';

%%
tor = rgb(:,:,1);
tog = rgb(:,:,2);
tob = rgb(:,:,3);

for i = 1: size(tor,1)
    for j = 1:size(tor,2)
        final = T * [tor(i,j);tog(i,j);tob(i,j);1];
        finalr(i,j) = final(1);
        finalg(i,j) = final(2);
        finalb(i,j) = final(3);
    end
end


%%
squarer = finalr(580:1564,3260:3917);
ch4r = mean(mean(squarer(824:956,512:640)));
checker4r = checker(4,1,1);
fr = checker4r/ch4r;
correctr = finalr * fr;

squareg = finalg(580:1564,3260:3917);
ch4g = mean(mean(squareg(824:956,512:640)));
checker4g = checker(4,1,2);
fg = checker4g/ch4g;
correctg = finalg * fg;

squareb = finalb(580:1564,3260:3917);
ch4b = mean(mean(squareb(824:956,512:640)));
checker4b = checker(4,1,3);
fb = checker4b/ch4b;
correctb = finalb * fb;

correct(:,:,1) = correctr;
correct(:,:,2) = correctg;
correct(:,:,3) = correctb;
%%
indexx = find(correct < 0);
correct(indexx) = 0;
hdrwrite(correct,'colorperfect.hdr');
hdrwrite(sqrt(correct/50),'color4.hdr');
%% Tone mapping method rgb separate
kkk = sqrt(correct/10);
K = 0.15;
B = 0.95;
N = size(kkk,1) * size(kkk,2);

%inter1 = kkk(:,:,1) + 1e-200; %r
%inter1 = kkk(:,:,2) + 1e-200; %g
inter1 = kkk(:,:,3) + 1e-200; %b

ImHDR = exp(1/N * sum(sum(log(inter1))));
I_ijHDR = K/ImHDR * kkk(:,:,1);
I_white = B * max(max(I_ijHDR));
I_ijTM = I_ijHDR .* (1+I_ijHDR/I_white/I_white)./(1+I_ijHDR);

%FINALIMG(:,:,1) = I_ijTM;
%FINALIMG(:,:,2) = I_ijTM;
FINALIMG(:,:,3) = I_ijTM;
 
%% Tone Mapping RGB as a whole
kkk = correct;
K = 3;
B = 3;
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
imwrite(FINALIMG,'rgbtonemapped.png');

%% Tone mapping with Y
kkk = correct;
K = 1;
B = 0.95;
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
imwrite(FINALIMG,'Y1tonemapped.png');




