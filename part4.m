cd ../data/exposure_stack
files = ["DSC_0299.jpg","DSC_0300.jpg",...
         "DSC_0301.jpg","DSC_0302.jpg","DSC_0303.jpg",...
         "DSC_0304.jpg","DSC_0305.jpg","DSC_0306.jpg",...
         "DSC_0307.jpg","DSC_0308.jpg","DSC_0309.jpg"];
%%
% crf = camresponse(files);  
expTimes = [1/512,1/256,1/128,1/64,1/32,1/16,1/8,1/4,1/2,1,2];
crf = camresponse(files,'ExposureTimes',expTimes);
range = 0:length(crf)-1;

figure,
hold on
plot(crf(:,1),range,'--r','LineWidth',2);
plot(crf(:,2),range,'-.g','LineWidth',2);
plot(crf(:,3),range,'-.b','LineWidth',2);
xlabel('Log-Exposure');
ylabel('Image Intensity');
title('Camera Response Function');
grid on
axis('tight')
legend('R-component','G-component','B-component','Location','southeast')


%%
img = imread('DSC_0302.jpg');
figure;
imshow(img)
title('Original file')
B = rgb2lin(img,'OutputType','double');
figure;
imshow(B)
title('Linearize with matlab function');
%%
img2 = imread('DSC_0306.jpg');
figure;
imshow(img2)
title('Original file')
C = rgb2lin(img2,'OutputType','double');
figure;
imshow(C)
title('Linearize with matlab function');

%%
cd ../../src
%% TOY
v1 = VideoReader('DSC_0284.MOV');
frame = read(v1,1);
figure(1)
imshow(frame)
title('Original');
L = rgb2lin(frame,'OutputType','double');
figure(2)
imshow(L)
title('Linearized');

v2 = VideoReader('DSC_0284.MOV');
frame2 = read(v2,400);
figure(3)
imshow(frame2)
title('Original with shadow');
LL = rgb2lin(frame2,'OutputType','double');
figure(4)
imshow(LL)
title('Linearized with shadow');
%%
v1 = VideoReader('DSC_0284.MOV');
num = 1;
for i = 300:5:450
    frame = read(v1,i);
    M = rgb2lin(frame,'OutputType','double');
    inter(:,:,:,num) = rgb2xyz(M,'ColorSpace','linear-rgb');
    Y(:,:,num) = inter(:,:,2);
    num = num+1;
    clear frame M
end

F = rgb2xyz(L,'ColorSpace','linear-rgb');
Y(:,:,num) = F(:,:,2);
inter(:,:,:,num) = F;

for r = 1:1080
    for c = 1:1920
        index_max(r,c) = min(find(Y(r,c,:) == max(Y(r,c,:))));
        index_min(r,c) = min(find(Y(r,c,:) == min(Y(r,c,:))));
    end
end

img_g = zeros(1080,1920,3);
for r = 1:1080
    for c = 1:1920
        img_g(r,c,:) = inter(r,c,:,index_min(r,c));
        img_d(r,c,:) = inter(r,c,:,index_max(r,c)) - inter(r,c,:,index_min(r,c));
    end
end
        
img_global = xyz2rgb(img_g,'ColorSpace','linear-rgb');
img_direct = xyz2rgb(img_d,'ColorSpace','linear-rgb');

figure(5)
imshow(img_global);

figure(6)
imshow(img_direct);




















%% BOTTLE
v1 = VideoReader('DSC_0282.MOV');
frame = read(v1,1);
figure(1)
imshow(frame)
title('Original');
L = rgb2lin(frame,'OutputType','double');
figure(2)
imshow(L)
title('Linearized');

v2 = VideoReader('DSC_0282.MOV');
frame2 = read(v2,400);
figure(3)
imshow(frame2)
title('Original with shadow');
LL = rgb2lin(frame2,'OutputType','double');
figure(4)
imshow(LL)
title('Linearized with shadow');
%%
num = 1;
for i = 400:5:550
    frame = read(v1,i);
    M = rgb2lin(frame,'OutputType','double');
    inter(:,:,:,num) = rgb2xyz(M,'ColorSpace','linear-rgb');
    Y(:,:,num) = inter(:,:,2);
    num = num+1;
    clear frame M
end

F = rgb2xyz(L,'ColorSpace','linear-rgb');
Y(:,:,num) = F(:,:,2);
inter(:,:,:,num) = F;

for r = 1:1080
    for c = 1:1920
        index_max(r,c) = min(find(Y(r,c,:) == max(Y(r,c,:))));
        index_min(r,c) = min(find(Y(r,c,:) == min(Y(r,c,:))));
    end
end

img_g = zeros(1080,1920,3);
for r = 1:1080
    for c = 1:1920
        img_g(r,c,:) = inter(r,c,:,index_min(r,c));
        img_d(r,c,:) = inter(r,c,:,index_max(r,c)) - inter(r,c,:,index_min(r,c));
    end
end
        
img_global = xyz2rgb(img_g,'ColorSpace','linear-rgb');
img_direct = xyz2rgb(img_d,'ColorSpace','linear-rgb');

figure(5)
imshow(img_global)
% title('Global');

figure(6)
imshow(img_direct)
% title('Direct');