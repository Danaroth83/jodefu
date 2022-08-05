function ATWT_Bspline(I_in, ratio, order)

u0=ones(1,ratio);
h=u0;
for ii=1:order-1
    h=conv(h,u0);
end
h=h/2^ratio;
Lh=length(h);

center_one=zeros(1,Lh);
center_one(ceil(Lh/2))=1;
g=center_one-h;

h_tilde=h;
g_tilde=center_one+h_tilde;

h=sqrt(2)*h;
g=sqrt(2)*g;
h_tilde=sqrt(2)*h_tilde;
g_tilde=sqrt(2)*g_tilde;


