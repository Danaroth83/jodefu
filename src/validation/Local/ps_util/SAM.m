function [angle_SAM,map] = SAM(msoriginal,msfused)

[M,N,~] = size(msfused);
prod_scal = dot(msoriginal,msfused,3); 
norm_orig = dot(msoriginal,msoriginal,3);
norm_fusa = dot(msfused,msfused,3);
prod_norm = sqrt(norm_orig.*norm_fusa);
prod_map = prod_norm;
prod_map(prod_map==0)=eps;
map = acos(prod_scal./prod_map);
prod_scal = reshape(prod_scal,M*N,1);
prod_norm = reshape(prod_norm, M*N,1);
z=find(prod_norm==0);
prod_scal(z)=[];prod_norm(z)=[];
angolo = sum(sum(acos(prod_scal./prod_norm)))/(size(prod_norm,1));
angle_SAM = real(angolo)*180/pi;

end