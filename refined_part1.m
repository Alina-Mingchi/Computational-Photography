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
Diff = Imax - Imin; 
DM = M;

% diff = Imax - Imin;
% [r c] = find(diff<100);

for i = 1:num_f
    DM(:,:,i) = I{i} - Ishadow; 
end

temp = DM;
for i = 1:num_f
    [rin cin] = find(Diff < 50);
    for j = 1: size(rin,1)
        temp(rin(j),cin(j),i) = 999;
    end
end


%% Debug
for num = 1:166
    test_up = temp(1:460,:,num);
    test_d = temp(461:end,:,num);
    test_up = imgaussfilt(test_up,2);
    test_d = imgaussfilt(test_d,2);

    [r1 c1] = find(abs(test_up) < 1 );
    [r2 c2] = find(abs(test_d) < 1 );
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
    pts(num,:) = [x1,y1,x2,y2,x3,y3,x4,y4];
%     figure(num);
%     imshow(M(:,:,num)/255)
%     hold on
%     plot(c1(left1),r1(left1),'g','LineWidth',3)
%     hold on
%     plot(c2(left2),r2(left2)+620,'g','LineWidth',3)
%     hold on
    
    clear c1 c2 r1 r2 left1 left2 FX grad_up grad_d test_up test_d
end
load('pts.mat');


%% shadow lines
%e.g. 136
load('../data/calib/Calib_Results.mat');
load('R_h.mat');
load('R_v.mat');
load('T_h.mat');
load('T_v.mat');
%p1 p2 in horizontal
% x = [670,736,0;752,746,0;338,730,0;148,774,0]';
for frame = 1:166
    if isempty(find(pts(frame,:) == 0))
    xx = pts(frame,:);
    x = [xx(1),xx(3),xx(5),xx(7);xx(2),xx(4),xx(6),xx(8)];
    r = pixel2ray(x,fc,cc,kc,alpha_c); % camera coordinate
    n1_c = R_v * [0;0;1];
    n2_c = R_v * [0;0;1];
    n3_c = R_h * [0;0;1];
    n4_c = R_h * [0;0;1];
    
    O1_c = R_v * [0;0;0] + T_v;
    O2_c = R_v * [0;0;0] + T_v;
    O3_c = R_h * [0;0;0] + T_h;
    O4_c = R_h * [0;0;0] + T_h;

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

    pts3d(frame,:) = [P1',P2',P3',P4'];
    end
end
save('pts3d.mat','pts3d')
%% shadow planes
for frame = 1:166
    if isempty(find(pts3d(frame,:) == 0))
        P1 = pts3d(frame,1:3);
        P2 = pts3d(frame,4:6);
        P3 = pts3d(frame,7:9);
        P4 = pts3d(frame,10:12);
        n_temp = cross((P2'-P1'),(P4'-P3'));
        n_hat = n_temp/norm(n_temp);
        shadow(frame,:) = [P3,n_hat'];
    end
end
save('shadow.mat','shadow')

%% Reconstruction
load('template.mat');
load('../data/calib/Calib_Results.mat');
template = template*160/255;
frame_num = template(326:640,314:808); %64-159
frog = M(326:640,314:808,:);
i=1;
for x = 326:640
    for y = 314:808
        v = pixel2ray([x;y],fc,cc,kc,alpha_c);
        if (floor(template(x,y))==template(x,y))
            if template(x,y) ~=0
                P(i,:) = shadow(template(x,y),1:3);
                n(i,:) = shadow(template(x,y),4:6);
                lambda(i) = dot(n(i,:)',P(i,:)')/dot(n(i,:)',v);
                Pttx(x-325,y-313) = lambda(i) * v(1);
                Ptty(x-325,y-313) = lambda(i) * v(2);
                Pttz(x-325,y-313) = lambda(i) * v(3);
            end
        end
        i = i+1;
    end
end
Ptt(:,:,1) = Pttx;
Ptt(:,:,2) = Ptty;
Ptt(:,:,3) = Pttz;

figure;
f000001 = im2double(imread('../data/frog/v1/000001.jpg'));
colormat = f000001(326:640,314:808,:);
ptCloud = pointCloud(Ptt,'Color', colormat);
% ptCloud.Color = colormat;
pcshow(ptCloud)


