% use the MATRIX from main.m
GG = {};
 for ii = 1:400
     for jj = 1:700
        MM = zeros(1,256);
        for count = 1:256
            aaa = mod(count,16);
            if aaa == 0
                aaa = 16;
            end
            temp = double(MAT{(count-aaa)/16+1,aaa});
            MM(count) = temp(ii,jj); 
        end
        GG = {GG,MM};
     end
 end
Mreal =  reshape(MM,[16,16]);



%%
for iii = 1:16
    for jjj = 1:16
        temp = double(MAT{iii,jjj});
        MMM(iii,jjj) = temp(100,300);
    end
end

figure()
imagesc(MMM)

%%
var = std(MMM).^2;
mini = min(var);
f = find(var == mini);






