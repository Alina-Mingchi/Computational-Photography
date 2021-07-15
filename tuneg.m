load('g.mat');
gg = g
for i = 1:255
    if (gg(i+1) < gg(i))
        gg(i+1) = gg(i);
    end
end
save('gg.mat','gg');
plot(1:256,gg)