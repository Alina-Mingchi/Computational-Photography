function [w numerator] = uniweight(image, zmin, zmax,t)
    [row col] = size(image);
    w = zeros(row, col);
    numerator = zeros(row, col);
    for i = 1:row
        for j = 1:col
            if (image(i,j) > zmin) && (image(i,j) < zmax)
                %uni weight
%                 w(i,j) = 1;
%                 numerator(i,j) = image(i,j)/t;
                
                %tent weight
                %w(i,j) = min(image(i,j),65535-image(i,j));
                %numerator(i,j) = w(i,j) .* image(i,j)/t;
                %tent jpeg
%                 w(i,j) = min(image(i,j),255-image(i,j));
%                 numerator(i,j) = w(i,j) .* image(i,j)/t;
                
                %gaussian weight
%                 w(i,j) = exp(-4*(image(i,j)-32768)*(image(i,j)-32768)/32768);
%                 numerator(i,j) = w(i,j) * image(i,j)/t;
                %gaussian jpeg
%                 w(i,j) = exp(-4*(image(i,j)-128)*(image(i,j)-128)/128);
%                 numerator(i,j) = w(i,j) * image(i,j)/t;
                
                %photon weight
%                 w(i,j) = t;
%                 numerator(i,j) = image(i,j);
%                 
                %log uni weight
%                 w(i,j) = 1;
%                 numerator(i,j) = log(image(i,j)/t);
                
                %log tent weight
                %w(i,j) = min(image(i,j),65535-image(i,j));
                %numerator(i,j) = w(i,j) .* log(image(i,j)/t);
                %tent jpeg
%                 w(i,j) = min(image(i,j),255-image(i,j));
%                 numerator(i,j) = w(i,j) .* log(image(i,j)/t);
                
                %log gaussian weight
%                 w(i,j) = exp(-4*(image(i,j)-32768)*(image(i,j)-32768)/32768);
%                 numerator(i,j) = w(i,j) * log(image(i,j)/t);
                %gaussian jpeg
%                 w(i,j) = exp(-4*(image(i,j)-128)*(image(i,j)-128)/128);
%                 numerator(i,j) = w(i,j) * log(image(i,j)/t);
                %photon weight
%                 w(i,j) = t;
%                 numerator(i,j) = w(i,j) * log(image(i,j)/t);
%part5
                w(i,j) = t^2/(33 * image(i,j)+1842);
                numerator(i,j) = w(i,j) * log(image(i,j)/t);
            end
        end
    end
end