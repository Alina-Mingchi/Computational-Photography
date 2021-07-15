
% f000001 = double(rgb2gray(imread('../data/frog/v1/000001.jpg')));

files=dir('../data/frog/v1/*.jpg');
for i=1:length(files);
    I{i}=double(rgb2gray(imread(files(i).name)));  
end

num_f = i; %166
for i = 1:num_f
    M(:,:,i) = I{i};
end
%%
Imax = max(M,[],3);
Imin = min(M,[],3);
Ishadow = (Imax + Imin)/2;
 
DM = M;

% diff = Imax - Imin;
% [r c] = find(diff<100);

for i = 1:num_f
    DM(:,:,i) = I{i} - Ishadow;    
end

%% plot on the original image

%% 136

%% patch1
temp1 = DM(:,:,136);
patch1 = temp1(140:180,740:820);
patch1 = imgaussfilt(patch1,2);
[r1 c1] = find(abs(patch1) < 2 );
figure;imshow(patch1/200);
hold on
plot(c1,r1,'*');
%%
FX = gradient(patch1);
for i = 1:size(r1,1)
    grad(i) = FX(r1(i),c1(i));
end

left1 = size(find(grad<0),2);
p1 = polyfit(c1(1:ceil(left1/2)),r1(1:ceil(left1/2)),1);

x1 = linspace(-50,150);
y1 = polyval(p1,x1);
figure;
imshow(patch1/100)
hold on
plot(c1(1:ceil(left1/2)),r1(1:ceil(left1/2)),'b*')
hold on
plot(x1,y1,'g')
hold on

right1 = size(find(grad<0),2);
pr1 = polyfit(c1(size(c1,1)-10:end),r1(size(c1,1)-10:end),1);

xr1 = linspace(0,200);
yr1 = polyval(pr1,xr1);
figure;
imshow(patch1/100)
hold on
plot(c1(size(c1,1)-10:end),r1(size(c1,1)-10:end),'r*')
hold on
plot(xr1,yr1,'b')
hold on

%%
patch2 = temp1(670:740,710:810);
patch2 = imgaussfilt(patch2,2);
[r1 c1] = find(abs(patch2) < 2 );
figure;imshow(patch2/200);
hold on
plot(c1,r1,'*');
%%
FX = gradient(patch2);
for i = 1:size(r1,1)
    grad(i) = FX(r1(i),c1(i));
end

left1 = size(find(grad<0),2);
pp1 = polyfit(c1(1:ceil(left1/2)),r1(1:ceil(left1/2)),1);

xx1 = linspace(-30,30);
yy1 = polyval(pp1,xx1);
figure;
imshow(patch2/100)
hold on
plot(c1(1:ceil(left1/2)),r1(1:ceil(left1/2)),'b*')
hold on
plot(xx1,yy1,'g')
hold on

right1 = size(find(grad<0),2);
ppr1 = polyfit(c1(size(c1,1)-10:end),r1(size(c1,1)-10:end),1);

xxr1 = linspace(-30,30);
yyr1 = polyval(ppr1,xxr1);
figure;
imshow(patch1/100)
hold on
plot(c1(size(c1,1)-10:end),r1(size(c1,1)-10:end),'r*')
hold on
plot(xxr1,yyr1,'b')
hold on

figure;
imshow(M(:,:,136)/max(max(M(:,:,136))));
hold on
plot(x1+740,y1+140,'g','LineWidth',3)
hold on
plot(xr1+740,yr1+140,'b','LineWidth',3)
hold on
plot(xx1+738,yy1+870,'g','LineWidth',3)
hold on
plot(xxr1+778,yyr1+1070,'b','LineWidth',3)
hold on


%% 90
% figure;
% imshow(M(:,:,90)/max(max(M(:,:,90))));
% hold on
% plot(x1+470,y1+170,'g','LineWidth',3)
% hold on
% plot(xr1+470,yr1+170,'b','LineWidth',3)
% hold on
% plot(x2+390,y2+1275,'g','LineWidth',3)
% hold on
% plot(xr2+495,yr2+4500,'b','LineWidth',3)
% hold on

%% 43

% figure;
% imshow(M(:,:,43)/max(max(M(:,:,43))));
% hold on
% plot(x1+250,y1+40,'g','LineWidth',3)
% hold on
% plot(xr1+250,yr1+40,'b','LineWidth',3)
% hold on


%% temporal
template = M(:,:,1);
indices = find((Imax - Imin)<50);
figure;
imshow(M(:,:,1)/255);
hold on;
for i = 20:160
    curr_img = DM(:,:,i);
    curr_img = imgaussfilt(curr_img,2);
    curr_img(indices) = 255;
    [r1 c1] = find(curr_img < 2 );
    plot(c1,r1,'*');
    hold on
    pix = find(curr_img < 2 );
    template(pix) = i*255/160;
end
figure;imagesc(template);colormap('jet')
save('template.mat','template');
