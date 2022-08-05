function [valori,cdf] = calc_cdf (v, num_interv)
v = reshape (v, 1, []);
inf= min(v);
sup= max(v);
lung_interv=(sup-inf)/num_interv;
cont=zeros(1,num_interv);
for i=1:length(v)
    for j=2:num_interv
        if v(i)>inf+(j-1)*lung_interv && v(i)<= inf+j*lung_interv
            cont(j)=cont(j)+1;
        end
    end
    if v(i) <= inf+lung_interv
        cont(1)=cont(1)+1;
    end
end
%asc=inf+lung_interv/2:lung_interv:sup-lung_interv/2;
cont = cont/length(v);% frequenza relativa
%bar(asc,cont);
cdf=[0 cumsum(cont)];
valori = inf:lung_interv:sup;
end