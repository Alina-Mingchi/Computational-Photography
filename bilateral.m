A = double(imread('../data/lamp/lamp_ambient.tif'))/255;
A_base = double(imread('../data/lamp/lamp_flash.tif'))/255;
figure(1)
imshow(A);

%For double checking the result
% J = imbilatfilt(A);
% figure(2)
% imshow(J);
% 
% figure(3)
% imshow(10 * abs(J-A));

%% bilateral filtering
[r c d] = size(A);
k_s = [];
sigma_s = 4; %f
sigma_r = 0.25; %g
hsize =4 * sigma_s + 1;

f = fspecial('gaussian',hsize,sigma_s);
g = fspecial('gaussian',hsize,sigma_r);
[ii jj] = meshgrid(1:hsize,1:hsize);
M = sqrt((ii-3).^2+(jj-3).^2);
fp_s = f .* M;

square = zeros(2*sigma_s);
rect1 = zeros(2*sigma_s,c);
rect2 = zeros(r,2*sigma_s);
Ip1 = [square,rect1,square;rect2,A(:,:,1),rect2;square,rect1,square];
Ip2 = [square,rect1,square;rect2,A(:,:,2),rect2;square,rect1,square];
Ip3 = [square,rect1,square;rect2,A(:,:,3),rect2;square,rect1,square];

%%
for j = 1:c
    for i = 1:r
        k_s(i,j,1) = sum(sum(fp_s*(g.*(Ip1(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,1)))));
        Js(i,j,1) = 1/k_s(i,j,1)*sum(sum(fp_s*(g.*(Ip1(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,1)))*Ip1(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end

%%
for j = 1:c
    for i = 1:r
        k_s(i,j,2) = sum(sum(fp_s*(g.*(Ip2(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,2)))));
        Js(i,j,2) = 1/k_s(i,j,2)*sum(sum(fp_s*(g.*(Ip2(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,2)))*Ip2(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end

%%
for j = 1:c
    for i = 1:r
        k_s(i,j,3) = sum(sum(fp_s*(g.*(Ip3(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,3)))));
        Js(i,j,3) = 1/k_s(i,j,3)*sum(sum(fp_s*(g.*(Ip3(i:i+4*sigma_s,j:j+4*sigma_s)-A(i,j,3)))*Ip3(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end

%%
Jss = Js/10;
TF = isnan(Jss);
Jss(TF) = A(TF);























%% Joint bilateral

sigma_s = 4; %f
sigma_r = 0.1; %g
hsize =4 * sigma_s + 1;

f = fspecial('gaussian',hsize,sigma_s);
g = fspecial('gaussian',hsize,sigma_r);
[ii jj] = meshgrid(1:hsize,1:hsize);
M = sqrt((ii-3).^2+(jj-3).^2);
fp_s = conv2(f,M);

square = zeros(2*sigma_s);
rect1 = zeros(2*sigma_s,c);
rect2 = zeros(r,2*sigma_s);
Ip1 = [square,rect1,square;rect2,A(:,:,1),rect2;square,rect1,square];
Ip2 = [square,rect1,square;rect2,A(:,:,2),rect2;square,rect1,square];
Ip3 = [square,rect1,square;rect2,A(:,:,3),rect2;square,rect1,square];

F1 = [square,rect1,square;rect2,A_base(:,:,1),rect2;square,rect1,square];
F2 = [square,rect1,square;rect2,A_base(:,:,2),rect2;square,rect1,square];
F3 = [square,rect1,square;rect2,A_base(:,:,3),rect2;square,rect1,square];

for j = 1:c
    for i = 1:r
        k_s(i,j,1) = sum(sum(fp_s*(conv2(g,(F1(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,1))))));
        Js(i,j,1) = 1/k_s(i,j,1)*sum(sum(fp_s*(conv2(g,(F1(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,1))))*Ip1(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end


for j = 1:c
    for i = 1:r
        k_s(i,j,2) = sum(sum(fp_s*(conv2(g,(F2(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,2))))));
        Js(i,j,2) = 1/k_s(i,j,2)*sum(sum(fp_s*(conv2(g,(F2(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,2))))*Ip2(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end


for j = 1:c
    for i = 1:r
        k_s(i,j,3) = sum(sum(fp_s*(conv2(g,(F3(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,3))))));
        Js(i,j,3) = 1/k_s(i,j,3)*sum(sum(fp_s*(conv2(g,(F3(i:i+4*sigma_s,j:j+4*sigma_s)-A_base(i,j,3))))*Ip3(i:i+4*sigma_s,j:j+4*sigma_s)));
    end
end










