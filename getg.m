%% read images and downsample them
%there used to be a subfolder within the data folder called exposure_stack, 
%where exist 16 RAW files, 16 JPEG files, and 16 TIFF files output from doing 
%dcraw in terminal, these are removed for the sake that the whole folder is too bulky
%That is why the path is shown as below

%% REad JPEG files
exo1 = imread('../data/exposure_stack/exposure1.jpg');
exo2 = imread('../data/exposure_stack/exposure2.jpg');
exo3 = imread('../data/exposure_stack/exposure3.jpg');
exo4 = imread('../data/exposure_stack/exposure4.jpg');
exo5 = imread('../data/exposure_stack/exposure5.jpg');
exo6 = imread('../data/exposure_stack/exposure6.jpg');
exo7 = imread('../data/exposure_stack/exposure7.jpg');
exo8 = imread('../data/exposure_stack/exposure8.jpg');
exo9 = imread('../data/exposure_stack/exposure9.jpg');
exo10 = imread('../data/exposure_stack/exposure10.jpg');
exo11 = imread('../data/exposure_stack/exposure11.jpg');
exo12 = imread('../data/exposure_stack/exposure12.jpg');
exo13 = imread('../data/exposure_stack/exposure13.jpg');
exo14 = imread('../data/exposure_stack/exposure14.jpg');
exo15 = imread('../data/exposure_stack/exposure15.jpg');
exo16 = imread('../data/exposure_stack/exposure16.jpg');
%%
N = 200;
esure1 = exo1(1:N:end,1:N:end);
esure2 = exo2(1:N:end,1:N:end);
esure3 = exo3(1:N:end,1:N:end);
esure4 = exo4(1:N:end,1:N:end);
esure5 = exo5(1:N:end,1:N:end);
esure6 = exo6(1:N:end,1:N:end);
esure7 = exo7(1:N:end,1:N:end);
esure8 = exo8(1:N:end,1:N:end);
esure9 = exo9(1:N:end,1:N:end);
esure10 = exo10(1:N:end,1:N:end);
esure11 = exo11(1:N:end,1:N:end);
esure12 = exo12(1:N:end,1:N:end);
esure13 = exo13(1:N:end,1:N:end);
esure14 = exo14(1:N:end,1:N:end);
esure15 = exo15(1:N:end,1:N:end);
esure16 = exo16(1:N:end,1:N:end);

%% linearize jpeg and normalize [0,1]
Z = [esure1(:),esure2(:),esure3(:),esure4(:),esure5(:),esure6(:),esure7(:),esure8(:),esure9(:),esure10(:),esure11(:),esure12(:),esure13(:),esure14(:),esure15(:),esure16(:)];
B = 1/2048 * [2^0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14,2^15];
lambda = 10; %try 0.1, 1, 10, 100
w = ones(257,1);

[g lE] = gsolve(Z,B,lambda,w);

figure
plot(1:256,g);
hold on
ylabel('g');
xlim([1 256]);
title('plot of g with lambda 10')

save('g.mat','g');