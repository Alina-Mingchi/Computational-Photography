%% part1

t = Tiff('../data/input_1.tif','r');
im1 = read(t);
t = Tiff('../data/input_2.tif','r');
im2 = read(t);
t = Tiff('../data/input_3.tif','r');
im3 = read(t);
t = Tiff('../data/input_4.tif','r');
im4 = read(t);
t = Tiff('../data/input_5.tif','r');
im5 = read(t);
t = Tiff('../data/input_6.tif','r');
im6 = read(t);
t = Tiff('../data/input_7.tif','r');
im7 = read(t);

lum1 = rgb2xyz(im1,'Colorspace','linear-rgb');
y1 = lum1(:,:,2);
lum2 = rgb2xyz(im2,'Colorspace','linear-rgb');
y2 = lum2(:,:,2);
lum3 = rgb2xyz(im3,'Colorspace','linear-rgb');
y3 = lum3(:,:,2);
lum4 = rgb2xyz(im4,'Colorspace','linear-rgb');
y4 = lum4(:,:,2);
lum5 = rgb2xyz(im5,'Colorspace','linear-rgb');
y5 = lum5(:,:,2);
lum6 = rgb2xyz(im6,'Colorspace','linear-rgb');
y6 = lum6(:,:,2);
lum7 = rgb2xyz(im7,'Colorspace','linear-rgb');
y7 = lum7(:,:,2);

[r c] = size(y7);

I = [y1(:)';y2(:)';y3(:)';y4(:)';y5(:)';y6(:)';y7(:)'];

load('../data/sources.mat');
L = S;

Be = pinv(L)*I;

for i = 1:size(Be,2)
    Ae(i) = norm(Be(:,i));
end

deno = repmat(Ae,3,1);
Ne = Be ./ deno;

A = reshape(Ae,[r,c]);
N(:,:,1) = reshape(Ne(1,:),[r,c]);
N(:,:,2) = reshape(Ne(2,:),[r,c]);
N(:,:,3) = reshape(Ne(3,:),[r,c]);

norm_N = (N+1)/2;
%%
figure(1);
imshow(A);

figure(2);
imshow(norm_N)

%% 
D = eye(3);
D(3,3) = -1;
Be = D * Be;

for i = 1:size(Be,2)
    Ae(i) = norm(Be(:,i));
end

deno = repmat(Ae,3,1);
Ne = Be ./ deno;

A = reshape(Ae,[r,c]);
N(:,:,1) = reshape(Ne(1,:),[r,c]);
N(:,:,2) = reshape(Ne(2,:),[r,c]);
N(:,:,3) = reshape(Ne(3,:),[r,c]);

fx = N(:,:,1) ./ N(:,:,3);
fy = N(:,:,2) ./ N(:,:,3);

Z = integrate_poisson(fx,fy);

maxi = max(max(Z));
mini = min(min(Z));

figure(12)
imshow(((Z-mini)/(maxi-mini)));


%%
figure;
s = surf(-Z);
axis image; axis off;
set(s,'facecolor',[.5 .5 .5],'edgecolor','none');
l = camlight;
rotate3d on





