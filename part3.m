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

%% 
load('pts.mat');

%% Shadow lines
load('../data/calib/Calib_Results.mat');
load('R_h.mat');
load('R_v.mat');
load('T_h.mat');
load('T_v.mat');


% ax + by + c = 0
for ii = 1:166
    parah(ii,1) = pts(ii,2) - pts(ii,4);
    parah(ii,2) = pts(ii,3) - pts(ii,1);
    parah(ii,3) = pts(ii,1) * pts(ii,4) - pts(ii,3) * pts(ii,2);
    parav(ii,1) = pts(ii,6) - pts(ii,8);
    parav(ii,2) = pts(ii,7) - pts(ii,5);
    parav(ii,3) = pts(ii,5) * pts(ii,8) - pts(ii,7) * pts(ii,6);
    
    if parah(ii,1) ~= 0
    
        syms x y
        eqns = [parah(ii,1)*x+parah(ii,2)*y+parah(ii,3) == 0, parav(ii,1)*x+parav(ii,2)*y+parav(ii,3) == 0];
        S = solve(eqns,[x y]);
        w(ii,1) = S.x;
        w(ii,2) = S.y;
        %w is shadow plane coordinate vector
        %if no intersection, w is (0,0)
    end
    
     o_v = R_v * [0;0;1] + T_v;
     o_h = R_h * [0;0;1] + T_h;

     l_v = [w(ii,:)';1] - o_v;
     l_h = [w(ii,:)';1] - o_h;
    
     A = [dot(l_h,l_h),-dot(l_h,l_v);-dot(l_h,l_v),dot(l_v,l_v)];
    b = [dot(l_h,(o_v-o_h));dot(l_v,(o_h-o_v))];
    result = A\b;
    alpha_h = result(1);
    alpha_v = result(2);
    w1(ii,:) = o_h + alpha_h * l_h;
    w2(ii,:) = o_v + alpha_v * l_v;
end
    



















