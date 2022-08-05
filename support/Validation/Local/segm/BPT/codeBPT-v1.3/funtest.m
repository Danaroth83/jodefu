H = [tree.nodeinfo];
H = [H.height];
m = max(H);

h = waitbar(0,'getting frames');
for i=1:m
    pt = pruneBPTheight(m-i+1);
    retrievesegmentation(pt,imrgb)
    F(i) = getframe;
    close
    waitbar(i/m)
end
close(h)    
%%
hf = figure;
set(hf,'Position',[100,100,size(F(1).cdata,2),size(F(1).cdata,1)])
movie(hf,F,1,4)
