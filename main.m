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

%% part2
[U,S,V] = svds(I);
U3 = U(:,1:3);
W3 = sqrt(S(1:3,1:3));
V3 = V(:,1:3)';

Le = U3 * W3;
Be = W3 * V3;

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
imshow(A*10);

figure(2);
imshow(norm_N)


%%

Q = [1,2,3;2,1,3;2,3,3];
Lq = Le * Q;
qqq = inv(Q);
Bq = qqq' * Be;
for j = 1:size(Bq,2)
    Aq(j) = norm(Bq(:,j));
end

deno = repmat(Aq,3,1);
Nq = Bq ./ deno;

AQ = reshape(Aq,[r,c]);
NQ(:,:,1) = reshape(Nq(1,:),[r,c]);
NQ(:,:,2) = reshape(Nq(2,:),[r,c]);
NQ(:,:,3) = reshape(Nq(3,:),[r,c]);

norm_NQ = (NQ+1)/2;
figure(3);
imshow(AQ*10);

figure(4);
imshow(norm_NQ)

%% part 3
l1 = [0.58,-0.58,-0.58];
l2 = [-0.58,-0.58,-0.58];
img_e1 = l1 * Be;
img_e2 = l2 * Be;
figure(5);
imshow(reshape(img_e1,[r,c])*10);
figure(6);
imshow(reshape(img_e2,[r,c])*10);

img_q1 = l1 * Bq;
img_q2 = l2 * Bq;
figure(7);
imshow(reshape(img_q1,[r,c])*10);
figure(8);
imshow(reshape(img_q2,[r,c])*10);
figure(9);
imshow((reshape(img_q2,[r,c])-min(img_q2))*10);

%% part 4
img(:,:,1) = reshape(Be(1,:),[r,c]);
img(:,:,2) = reshape(Be(2,:),[r,c]);
img(:,:,3) = reshape(Be(3,:),[r,c]);

newBe = imgaussfilt(img,20);
[FX,FY] = gradient(newBe);


temp1 = newBe(:,:,1);
temp2 = newBe(:,:,2);
temp3 = newBe(:,:,3);
tempp1= temp1(:);
tempp2= temp2(:);
tempp3= temp3(:);
newBe = [tempp1';tempp2';tempp3'];

temp1 = FX(:,:,1);
temp2 = FX(:,:,2);
temp3 = FX(:,:,3);
tempp1= temp1(:);
tempp2= temp2(:);
tempp3= temp3(:);
FXX = [tempp1';tempp2';tempp3'];

temp1 = FY(:,:,1);
temp2 = FY(:,:,2);
temp3 = FY(:,:,3);
tempp1= temp1(:);
tempp2= temp2(:);
tempp3= temp3(:);
FYY = [tempp1';tempp2';tempp3'];

%%
A1 = newBe(1,:).*FXX(2,:) - newBe(2,:).*FXX(1,:);
A2 = newBe(1,:).*FXX(3,:) - newBe(3,:).*FXX(1,:);
A3 = newBe(2,:).*FXX(3,:) - newBe(3,:).*FXX(2,:);
A4 = -newBe(1,:).*FYY(2,:) + newBe(2,:).*FYY(1,:);
A5 = -newBe(1,:).*FYY(3,:) + newBe(3,:).*FYY(1,:);
A6 = -newBe(2,:).*FYY(3,:) + newBe(3,:).*FYY(2,:);

AAA = [A1;A2;A3;A4;A5;A6]';
[uu ss vv] = svds(AAA);
x = vv(:,6);
D = [-x(3) x(6) 1;x(2) -x(5) 0;-x(1) x(4) 0];
Btt = inv(D) * Be;
% alpha = 0.5;
% beta = 0.5;
alpha = 0;
beta = 0;
G = [1 0 0;0 1 0;alpha beta -1];
Btt = inv(G') * Btt;

albedoes = sqrt(sum(Btt .^ 2, 1));

normals = Btt ./ repmat(albedoes, [3, 1]);

NNN = normals;

ATT = reshape(albedoes,[r,c]);
NTT(:,:,1) = reshape(NNN(1,:),[r,c]);
NTT(:,:,2) = reshape(NNN(2,:),[r,c]);
NTT(:,:,3) = reshape(NNN(3,:),[r,c]);

%%

norm_NTT = (NTT+1)/2;
figure(10);
imshow(ATT*10);

figure(11);
imshow(norm_NTT);

%% part 5

%[fx,fy] = gradient(inter);
fx = NTT(:,:,1) ./ NTT(:,:,3);
fy = NTT(:,:,2) ./ NTT(:,:,3);

Z = integrate_poisson(fx,fy);
% Z = integrate_frankot(fx,fy);

% sigma = 10
% Z(Z>100) =100;
% Z(Z<-100)=-100;

figure(12)
imshow((Z+100)/200)

% %sigma = 1
% % Z(Z>200) =200;
% % Z(Z<-200)=-200;
% % 
% % figure(12)
% % imshow((Z+200)/400)


%%
figure(20);
s = surf(-Z);
axis image; axis off;
set(s,'facecolor',[.5 .5 .5],'edgecolor','none');
l = camlight;
rotate3d on










