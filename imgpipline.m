%% Matlab initials
image = imread('banana_slug.tiff');
[r,c] = size(image);
className = class(image);
whos image
img_double = double(image);
whos img_double

%% Linearization
maximum = max(max(img_double));
minimum = min(min(img_double));

trans1 = 2047 * ones(r,c);
trans2 = (15000 - 2047) * ones(r,c);
img_temp = (img_double - trans1) ./ trans2;
img_linear = img_temp;

for i = 1:r
    for j = 1:c
        if img_temp(i,j) < 0
            img_linear(i,j) = 0;
        elseif img_temp(i,j) > 1
            img_linear(i,j) = 1;
        end
    end
end

% figure;
% imshow(min(1,img_linear * 5));
    
%% Identify Bayer pattern
im1 = img_linear(1:2:end,1:2:end);
av1 = mean(mean(im1));
figure(1);
imshow(min(1,im1 * 5));
im2 = img_linear(1:2:end,2:2:end);
av2 = mean(mean(im2));
figure(2);
imshow(min(1,im2 * 5));
im3 = img_linear(2:2:end,1:2:end);
av3 = mean(mean(im3));
figure(3);
imshow(min(1,im3 * 5));
im4 = img_linear(2:2:end,2:2:end);
av4 = mean(mean(im4));
figure(4);
imshow(min(1,im4 * 5));    
    
%% White balancing    
% grey world assumption
avr = av1;
avg = (av2 + av3)/2;
avb = av4;

r_p = avg/avr * ones(r/2,c/2) .* im1;
g_p1 = av2/avg * ones(r/2,c/2) .* im2;
g_p2 = av3/avg * ones(r/2,c/2) .* im3;
b_p = avb/avr * ones(r/2,c/2) .* im4;

%white world assumption
r_max = max(max(im1));
g1_max = max(max(im2));
g2_max = max(max(im3));
g_max = max(g1_max, g2_max);
b_max = max(max(im4));

r_pp = g_max/r_max* ones(r/2,c/2) .* im1;
g_pp1 = g_max/g1_max * ones(r/2,c/2) .* im2;
g_pp2 = g_max/g2_max * ones(r/2,c/2) .* im3;
b_pp = g_max/b_max * ones(r/2,c/2) .* im4;

%% Demosaicing
rrow = size(im1,1);
rows = 1:rrow;
rcol = size(im1,2);
columns = 1:rcol;
[X Y]=meshgrid(columns,rows);
[Xq Yq]=meshgrid(columns,rows);
V = im1;
Vr = interp2(X,Y,V,Xq,Yq,'cubic');

brow = size(im4,1);
rows = 1:brow;
bcol = size(im4,2);
columns = 1:bcol;
[X Y]=meshgrid(columns,rows);
[Xq Yq]=meshgrid(columns,rows);
V = im4;
Vb = interp2(X,Y,V,Xq,Yq,'cubic');

grow = size(im2,1);
rows = 1:grow;
gcol = size(im2,2);
columns = 1:gcol;
[X Y]=meshgrid(columns,rows);
[Xq Yq]=meshgrid(columns,rows);
V = (im2+im3)/2;
Vg = interp2(X,Y,V,Xq,Yq,'cubic');

%When generating gray world assumpltion, uncomment the following lines
%bilinear interpolation by hand
% %gray
% A = r_p;
% B = zeros(r/2,1);
% B = [B,A];
% AB = A + B(:,2:c/2+1);
% temp = zeros(1,c/2);
% temp = [temp;AB];
% red = (AB + temp(2:r/2+1,:))/4;
% 
% C = b_p;
% D = zeros(r/2,1);
% D = [D,C];
% CD = C + D(:,2:c/2+1);
% tempp = zeros(1,c/2);
% tempp = [tempp;CD];
% blue = (CD + tempp(2:r/2+1,:))/4;
% 
% G1 = g_p1;
% E = zeros(1,c/2);
% E = [E;G1];
% greenn= G1 + E(2:r/2+1,:);
% G2 = g_p2;
% F = zeros(r/2,1);
% F = [F,G2];
% greennn = G2 + F(:,2:c/2+1);
% green = (greenn + greennn)/4;

% white
A = r_pp;
B = zeros(r/2,1);
B = [B,A];
AB = A + B(:,2:c/2+1);
temp = zeros(1,c/2);
temp = [temp;AB];
red = (AB + temp(2:r/2+1,:))/4;

C = b_pp;
D = zeros(r/2,1);
D = [D,C];
CD = C + D(:,2:c/2+1);
tempp = zeros(1,c/2);
tempp = [tempp;CD];
blue = (CD + tempp(2:r/2+1,:))/4;

G1 = g_pp1;
E = zeros(1,c/2);
E = [E;G1];
greenn= G1 + E(2:r/2+1,:);
G2 = g_pp2;
F = zeros(r/2,1);
F = [F,G2];
greennn = G2 + F(:,2:c/2+1);
green = (greenn + greennn)/4;
%% Brightness adjustment & gamma correction 
intermediate = zeros(r/2,c/2,3);
intermediate(:,:,1) = red;
intermediate(:,:,2) = green;
intermediate(:,:,3) = blue;
I = rgb2gray(intermediate);
max_gray = max(max(I));
factor = 1/max_gray*2;
final = intermediate * factor;

for i = 1:r/2
    for j = 1:c/2
        for k = 1:3
            if final(i,j,k) <= 0.0031308
                final(i,j,k) = 12.92 * final(i,j,k);
            else
                final(i,j,k) = (1+0.055)*final(i,j,k)^(1/2.4)-0.055;
            end
        end
    end
end

%% Compression
% imwrite(final,'greyassumption.png');
imwrite(final,'brightwhiteassumption.png');
% imwrite(final,'whiteassumption.jpeg','jpeg','quality',95);
% imwrite(final,'whiteassumption10.jpeg','jpeg','quality',10);
% imwrite(final,'whiteassumption20.jpeg','jpeg','quality',20);
% imwrite(final,'whiteassumption30.jpeg','jpeg','quality',30);
% imwrite(final,'whiteassumption40.jpeg','jpeg','quality',40);
% imwrite(final,'whiteassumption50.jpeg','jpeg','quality',50);
% imwrite(final,'whiteassumption60.jpeg','jpeg','quality',60);
% imwrite(final,'whiteassumption70.jpeg','jpeg','quality',70);
% imwrite(final,'whiteassumption80.jpeg','jpeg','quality',80);
%% Manual white balancing
im1 = img_linear(1:2:end,1:2:end);
av1 = mean(mean(im1));
figure(1);
imshow(min(1,im1 * 5));
impixelinfo

%% small patch
smallr = im1(170:190,1030:1050);
smallg = (im2(170:190,1030:1050) + im3(170:190,1030:1050))/2;
smallb = im4(170:190,1030:1050);

smallraverage = mean(mean(smallr));
smallgaverage = mean(mean(smallg));
smallbaverage = mean(mean(smallb));

fac_r = smallgaverage/smallraverage;
fac_b = smallgaverage/smallbaverage;

r_pp = fac_r * im1;
b_pp = fac_b * im4;
g_pp1 = im2;
g_pp2 = im3;

A = r_pp;
B = zeros(r/2,1);
B = [B,A];
AB = A + B(:,2:c/2+1);
temp = zeros(1,c/2);
temp = [temp;AB];
red = (AB + temp(2:r/2+1,:))/4;

C = b_pp;
D = zeros(r/2,1);
D = [D,C];
CD = C + D(:,2:c/2+1);
tempp = zeros(1,c/2);
tempp = [tempp;CD];
blue = (CD + tempp(2:r/2+1,:))/4;

G1 = g_pp1;
E = zeros(1,c/2);
E = [E;G1];
greenn= G1 + E(2:r/2+1,:);
G2 = g_pp2;
F = zeros(r/2,1);
F = [F,G2];
greennn = G2 + F(:,2:c/2+1);
green = (greenn + greennn)/4;

intermediate = zeros(r/2,c/2,3);
intermediate(:,:,1) = red;
intermediate(:,:,2) = green;
intermediate(:,:,3) = blue;
I = rgb2gray(intermediate);
max_gray = max(max(I));
factor = 1/max_gray*2;
final = intermediate * factor;

for i = 1:r/2
    for j = 1:c/2
        for k = 1:3
            if final(i,j,k) <= 0.0031308
                final(i,j,k) = 12.92 * final(i,j,k);
            else
                final(i,j,k) = (1+0.055)*final(i,j,k)^(1/2.4)-0.055;
            end
        end
    end
end

imwrite(final,'manualbrightwhiteassumption.png');
%% large patch
smallr = im1(940:980,1320:1360);
smallg = (im2(940:980,1320:1360) + im3(940:980,1320:1360))/2;
smallb = im4(940:980,1320:1360);

smallraverage = mean(mean(smallr));
smallgaverage = mean(mean(smallg));
smallbaverage = mean(mean(smallb));

fac_r = smallgaverage/smallraverage;
fac_b = smallgaverage/smallbaverage;

r_pp = fac_r * im1;
b_pp = fac_b * im4;
g_pp1 = im2;
g_pp2 = im3;

A = r_pp;
B = zeros(r/2,1);
B = [B,A];
AB = A + B(:,2:c/2+1);
temp = zeros(1,c/2);
temp = [temp;AB];
red = (AB + temp(2:r/2+1,:))/4;

C = b_pp;
D = zeros(r/2,1);
D = [D,C];
CD = C + D(:,2:c/2+1);
tempp = zeros(1,c/2);
tempp = [tempp;CD];
blue = (CD + tempp(2:r/2+1,:))/4;

G1 = g_pp1;
E = zeros(1,c/2);
E = [E;G1];
greenn= G1 + E(2:r/2+1,:);
G2 = g_pp2;
F = zeros(r/2,1);
F = [F,G2];
greennn = G2 + F(:,2:c/2+1);
green = (greenn + greennn)/4;

intermediate = zeros(r/2,c/2,3);
intermediate(:,:,1) = red;
intermediate(:,:,2) = green;
intermediate(:,:,3) = blue;
I = rgb2gray(intermediate);
max_gray = max(max(I));
factor = 1/max_gray*2;
final = intermediate * factor;

for i = 1:r/2
    for j = 1:c/2
        for k = 1:3
            if final(i,j,k) <= 0.0031308
                final(i,j,k) = 12.92 * final(i,j,k);
            else
                final(i,j,k) = (1+0.055)*final(i,j,k)^(1/2.4)-0.055;
            end
        end
    end
end

imwrite(final,'manual2brightwhiteassumption.png');


