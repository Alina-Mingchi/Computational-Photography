%% read images and downsample them
%there used to be a subfolder within the data folder called exposure_stack, 
%where exist 16 RAW files, 16 JPEG files, and 16 TIFF files output from doing 
%dcraw in terminal, these are removed for the sake that the whole folder is too bulky
%That is why the path is shown as below

%% REad JPEG files
expo1 = double(imread('../data/exposure_stack/exposure1.jpg'));
expo2 = double(imread('../data/exposure_stack/exposure2.jpg'));
expo3 = double(imread('../data/exposure_stack/exposure3.jpg'));
expo4 = double(imread('../data/exposure_stack/exposure4.jpg'));
expo5 = double(imread('../data/exposure_stack/exposure5.jpg'));
expo6 = double(imread('../data/exposure_stack/exposure6.jpg'));
expo7 = double(imread('../data/exposure_stack/exposure7.jpg'));
expo8 = double(imread('../data/exposure_stack/exposure8.jpg'));
expo9 = double(imread('../data/exposure_stack/exposure9.jpg'));
expo10 = double(imread('../data/exposure_stack/exposure10.jpg'));
expo11 = double(imread('../data/exposure_stack/exposure11.jpg'));
expo12 = double(imread('../data/exposure_stack/exposure12.jpg'));
expo13 = double(imread('../data/exposure_stack/exposure13.jpg'));
expo14 = double(imread('../data/exposure_stack/exposure14.jpg'));
expo15 = double(imread('../data/exposure_stack/exposure15.jpg'));
expo16 = double(imread('../data/exposure_stack/exposure16.jpg'));

%%
esure1 = [expo1(:,:,1),expo1(:,:,2),expo1(:,:,3)];
esure2 = [expo2(:,:,1),expo2(:,:,2),expo2(:,:,3)];
esure3 = [expo3(:,:,1),expo3(:,:,2),expo3(:,:,3)];
esure4 = [expo4(:,:,1),expo4(:,:,2),expo4(:,:,3)];
esure5 = [expo5(:,:,1),expo5(:,:,2),expo5(:,:,3)];
esure6 = [expo6(:,:,1),expo6(:,:,2),expo6(:,:,3)];
esure7 = [expo7(:,:,1),expo7(:,:,2),expo7(:,:,3)];
esure8 = [expo8(:,:,1),expo8(:,:,2),expo8(:,:,3)];
esure9 = [expo9(:,:,1),expo9(:,:,2),expo9(:,:,3)];
esure10 = [expo10(:,:,1),expo10(:,:,2),expo10(:,:,3)];
esure11 = [expo11(:,:,1),expo11(:,:,2),expo11(:,:,3)];
esure12 = [expo12(:,:,1),expo12(:,:,2),expo12(:,:,3)];
esure13 = [expo13(:,:,1),expo13(:,:,2),expo13(:,:,3)];
esure14 = [expo14(:,:,1),expo14(:,:,2),expo14(:,:,3)];
esure15 = [expo15(:,:,1),expo15(:,:,2),expo15(:,:,3)];
esure16 = [expo16(:,:,1),expo16(:,:,2),expo16(:,:,3)];


%%
load('g.mat');
[r c] = size(esure1);
exposure1 = zeros(r,c);
exposure2 = zeros(r,c);
exposure3 = zeros(r,c);
exposure4 = zeros(r,c);
exposure5 = zeros(r,c);
exposure6 = zeros(r,c);
exposure7 = zeros(r,c);
exposure8 = zeros(r,c);
exposure9 = zeros(r,c);
exposure10 = zeros(r,c);
exposure11 = zeros(r,c);
exposure12 = zeros(r,c);
exposure13 = zeros(r,c);
exposure14 = zeros(r,c);
exposure15 = zeros(r,c);
exposure16 = zeros(r,c);
for i = 1:256
    index1 = find(esure1 == (i-1));
    exposure1(index1) = g(i);
    index2 = find(esure2 == (i-1));
    exposure2(index2) = g(i);
    index3 = find(esure3 == (i-1));
    exposure3(index3) = g(i);
    index4 = find(esure4 == (i-1));
    exposure4(index4) = g(i);
end

%%
for i = 1:256
    index16 = find(esure16 == (i-1));
    exposure16(index16) = g(i);
end
save('l16.mat','exposure16')

%% JPEG linear uniform
zmin = 0.01 * 255;
zmax = 0.99 * 255;
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
%% Linear Uniform
I_j_l_u = numerator ./ denominator;
save('j_li_u.mat','I_j_l_u');
%% Linear Tent
I_j_l_t = numerator ./ denominator;
save('j_li_t.mat','I_j_l_t');
%% Linear gaussian
I_j_l_g = numerator ./ denominator;
save('j_li_g.mat','I_j_l_g');
%% Linear photon
I_j_l_p = numerator ./ denominator;
save('j_li_p.mat','I_j_l_p');
%% Log uni
I_j_lo_u = exp(numerator ./ denominator);
save('j_lo_u.mat','I_j_lo_u');
%% Log tent
I_j_lo_t = exp(numerator ./ denominator);
save('j_lo_t.mat','I_j_lo_t');
%% Log gaussian
I_j_lo_g = exp(numerator ./ denominator);
save('j_lo_g.mat','I_j_lo_g');
%% Log photon
I_j_lo_p = exp(numerator ./ denominator);
save('j_lo_p.mat','I_j_lo_p');

%% Evaluation
% rgb(:,:,1) = I_j_l_u(:,1:6000);
% rgb(:,:,2) = I_j_l_u(:,6001:12000);
% rgb(:,:,3) = I_j_l_u(:,12001:end);
% 
% rgb(:,:,1) = I_j_l_t(:,1:6016);
% rgb(:,:,2) = I_j_l_t(:,6017:12032);
% rgb(:,:,3) = I_j_l_t(:,12033:end);

% rgb(:,:,1) = I_j_l_g(:,1:6016);
% rgb(:,:,2) = I_j_l_g(:,6017:12032);
% rgb(:,:,3) = I_j_l_g(:,12033:end);
%  
% rgb(:,:,1) = I_j_l_p(:,1:6016);
% rgb(:,:,2) = I_j_l_p(:,6017:12032);
% rgb(:,:,3) = I_j_l_p(:,12033:end);

% rgb(:,:,1) = I_j_lo_u(:,1:6016);
% rgb(:,:,2) = I_j_lo_u(:,6017:12032);
% rgb(:,:,3) = I_j_lo_u(:,12033:end);
% 
% rgb(:,:,1) = I_j_lo_t(:,1:6016);
% rgb(:,:,2) = I_j_lo_t(:,6017:12032);
% rgb(:,:,3) = I_j_lo_t(:,12033:end);

% rgb(:,:,1) = I_j_lo_g(:,1:6016);
% rgb(:,:,2) = I_j_lo_g(:,6017:12032);
% rgb(:,:,3) = I_j_lo_g(:,12033:end);
% % 
rgb(:,:,1) = I_j_lo_g(:,1:6000);
rgb(:,:,2) = I_j_lo_g(:,6001:12000);
rgb(:,:,3) = I_j_lo_g(:,12001:end);

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
title('plot of log regression of JPEG photon weights')
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









