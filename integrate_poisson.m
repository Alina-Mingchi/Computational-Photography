function Z = integrate_poisson(zx, zy)
% Least squares solution
% Poisson Reconstruction Using Neumann boundary conditions
% Input gx and gy
% gradients at boundary are assumed to be zero
% Output : reconstruction

[h, w] = size(zx);

zx(:,end) = 0;
zy(end,:) = 0;

% pad
zx = padarray(zx,[1 1],0,'both');
zy = padarray(zy,[1 1],0,'both');

gxx = zeros(size(zx)); 
gyy = zeros(size(zx)); 

j = 1:h+1;
k = 1:w+1;

% Laplacian
gyy(j+1,k) = zy(j+1,k) - zy(j,k);
gxx(j,k+1) = zx(j,k+1) - zx(j,k);
f = gxx + gyy;
f = f(2:end-1,2:end-1);
clear j k gxx gyy gyyd gxxd gx gy

%compute eigen values
[x,y] = meshgrid(0:w-1,0:h-1);
denom = (2*cos(pi*x/(w))-2) + (2*cos(pi*y/(h)) - 2);
clear x y
	
%compute cosine transform
fcos = dct2(f);
clear f
		
%divide. 1st element of denom will be zero and 1st element of fcos and
%after division should also be zero; so divided rest of the elements
fcos(2:end) = fcos(2:end)./denom(2:end);
clear denom

% compute inverse dct2
Z = idct2(fcos);

end
