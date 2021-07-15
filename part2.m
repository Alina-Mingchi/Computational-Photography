v = VideoReader('DSC_0273.MOV');
frame = read(v,1);
figure
imshow(frame)
img = rgb2gray(im2double(frame));
figure;imshow(img);
%%
num = 1;
for i = 500:5:700
    frame = read(v,i);
%     figure
%     imshow(frame(:,:,:,i))
    M(:,:,num) = rgb2gray(im2double(frame));
    num = num+1;
end

%%
Imax = max(M,[],3);
Imin = min(M,[],3);
Ishadow = (Imax + Imin)/2;
 
DM = M;

% diff = Imax - Imin;
% [r c] = find(diff<100);

for i = 1:41
    DM(:,:,i) = M(:,:,i) - Ishadow;    
end


%% temporal
template = img;
indices = find((Imax - Imin)<0.1);
figure;
imshow(template);
hold on;
i = 15;
for i = 1:41
    curr_img = DM(:,:,i);
    curr_img = imgaussfilt(curr_img,2);
    curr_img(indices) = 1;
    [r1 c1] = find(curr_img < -0.1 );
    plot(c1,r1,'*');
    hold on
    pix = find(curr_img < -0.1 );
    template(pix) = i;
end
% figure;imagesc(template*6);colormap('jet')


%%
%% Debug
temp = M;
for num = 1:41
    test_up = temp(1:460,:,num);
    test_d = temp(600:end,:,num);
    test_up = imgaussfilt(test_up,2);
    test_d = imgaussfilt(test_d,2);

    [r1 c1] = find(abs(test_up) < 1 );
    [r2 c2] = find(abs(test_d) < -0.1 );
% figure;imshow(M(:,:,100)/255);
% hold on
% plot(c1,r1,'.');
% hold on
% plot(c2,r2,'.');
    grad_up = [];
    grad_d = [];
    FX = gradient(test_up);
    for i = 1:size(r1,1)
        grad_up(i) = FX(r1(i),c1(i));
    end

    clear FX;

    FX = gradient(test_d);
    for i = 1:size(r2,1)
        grad_d(i) = FX(r2(i),c2(i));
    end

    left1 = find(grad_up<0);%v
    left2 = find(grad_d<0);%h
    
    size1 = size(left1,2);
    size2 = size(left2,2);
    
    if size2 == 0
        x1 = 0;
        y1 = 0;
        x2 = 0;
        y2 = 0;
    else
        x1 = c2(left2(size2));
        y1 = r2(left2(size2));
        x2 = c2(left2(ceil(size2/4)));
        y2 = r2(left2(ceil(size2/4)));
    end
    if size1 ==0
        x3 = 0;
        y3 = 0;
        x4 = 0;
        y4 = 0;
    else
        x3 = c1(left1(size1));
        y3 = r1(left1(size1));
        x4 = c1(left1(ceil(size1/4)));
        y4 = r1(left1(ceil(size1/4)));
    end
    ptts(num,:) = [x1,y1,x2,y2];
    
    clear c1 c2 r1 r2 left1 left2 FX grad_up grad_d test_up test_d
end
 load('ptts.mat');


%% shadow lines
load('../data/cal/Calib_Results.mat');
load('Rh.mat');
load('Rv.mat');
load('Th.mat');
load('Tv.mat');
%p1 p2 in horizontal
for frame = 1:41
    if isempty(find(ptts(frame,:) == 0))
    xx = ptts(frame,:);
    x = [xx(1),xx(3),xx(5),xx(7);xx(2),xx(4),xx(6),xx(8)];
    r = pixel2ray(x,fc,cc,kc,alpha_c); % camera coordinate
    n3_c = Rv * [0;0;1];
    n4_c = Rv * [0;0;1];
    n1_c = Rh * [0;0;1];
    n2_c = Rh * [0;0;1];
    
    O3_c = Rv * [0;0;0] + Tv;
    O4_c = Rv * [0;0;0] + Tv;
    O1_c = Rh * [0;0;0] + Th;
    O2_c = Rh * [0;0;0] + Th;

    A1 = dot(n1_c,O1_c);
    B1 = dot(n1_c,r(:,1));
    lambda1 = A1/B1;
    P1 = lambda1 * r(:,1);

    A2 = dot(n2_c,O2_c);
    B2 = dot(n2_c,r(:,2));
    lambda2 = A2/B2;
    P2 = lambda2 * r(:,2);

    A3 = dot(n3_c,O3_c);
    B3 = dot(n3_c,r(:,3));
    lambda3 = A3/B3;
    P3 = lambda3 * r(:,3);

    A4 = dot(n4_c,O4_c);
    B4 = dot(n4_c,r(:,4));
    lambda4 = A4/B4;
    P4 = lambda4 * r(:,4);

    ptts3d(frame,:) = [P1',P2',P3',P4'];
    end
end
save('ptts3d.mat','ptts3d')
%% shadow planes
for frame = 1:41
    if isempty(find(ptts3d(frame,:) == 0))
        P1 = ptts3d(frame,1:3);
        P2 = ptts3d(frame,4:6);
        P3 = ptts3d(frame,7:9);
        P4 = ptts3d(frame,10:12);
        n_temp = cross((P2'-P1'),(P4'-P3'));
%         n_hat = n_temp/norm(n_temp);
        sshadow(frame,:) = [P3,n_temp'];
    end
end
save('sshadow.mat','sshadow')

%% Reconstruction

load('../data/cal/Calib_Results.mat');
template = template;
frame_num = [];
frame_num = template(566:1060,370:1086); 
% box = M(566:1060,370,1086,:);
k = 1;
figure; hold on
% for x = 866:1060
%     for y = 570:1086
for x = 566:1060
    for y = 370:1086
        v = pixel2ray([x;y],fc,cc,kc,alpha_c);
        if (floor(template(x,y))==template(x,y))
            if template(x,y) ~=0
                P(k,:) = sshadow(template(x,y),1:3);
                n(k,:) = sshadow(template(x,y),4:6);
                lambda(k) = dot(n(k,:)',P(k,:)')/dot(n(k,:)',v);
                Pttx(x-565,y-369) = lambda(k) * v(1);
                Ptty(x-565,y-369) = lambda(k) * v(2);
                Pttz(x-565,y-369) = lambda(k) * v(3);
                
%                 plot3(lambda(k) * v(1),lambda(k) * v(2),lambda(k) * v(3),'.')
            end
        end
        k = k+1;
    end
end
Pttt(:,:,1) = Pttx;
Pttt(:,:,2) = Ptty;
Pttt(:,:,3) = Pttz;

figure(10);
vv = VideoReader('DSC_0273.MOV');
frame1 = read(vv,1);
img = im2double(frame1);
colormat =img(566:1060,370:1086,:);
ptCloud = pointCloud(Pttt,'Color', colormat);
pcshow(ptCloud)






