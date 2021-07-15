function Z = integrate_frankot(zx, zy)
% Integration using the Frankot-Chellappa algorithm.
% Input zx and zy
% Output : reconstruction

[h, w]=size(zx);

% complain if P or Q are too big
if (h>512) | (w>512)
  error('Input array too big.  Choose a smaller window.');
end

% pad the input array to 512x512
nrows=2^9; ncols=2^9;
   
% compute Fourier coefficients
if isempty(nrows)
   Zx=fft2(zx);
   Zy=fft2(zy);
   h2=h;w2=w;
else
   Zx=fft2(zx,nrows,ncols);
   Zy=fft2(zy,nrows,ncols);
   h2=nrows; w2=ncols;
end
Zx=Zx(:); Zy=Zy(:);

% compute repeated frequency vectors (See Chellapa paper)
Wy=repmat(2*pi/h2*[0:(h2/2), (-(h2/2)+1):(-1)]',w2,1);
Wx=kron(2*pi/w2*[0:(w2/2), (-(w2/2)+1):(-1)]',ones(h2,1));

% compute transform of least squares closest integrable surface
%    remove first column because it's all zeros (then add C(0)=0)
C=(-j*Wx(2:end).*Zx(2:end)-j*Wy(2:end).*Zy(2:end))./...
  (Wx(2:end).^2+Wy(2:end).^2);

% set DC component of C
C=[0;C];

% invert transform to get depth of integrable surface
Z=real(ifft2(reshape(C,h2,w2)));

% crop output if there was padding         
if ~isempty(nrows)
   Z=Z(1:h,1:w);
end   








