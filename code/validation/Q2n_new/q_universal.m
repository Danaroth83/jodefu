function [qs,qi,qj,qk]=q_universal(a1,b1,c1,d1,a2,b2,c2,d2,size1,size2);
% q2=a2+ib2+jc2+kd2
% q1=a1+ib1+jc1+kd1
% (size1,size2) is the dimension of the block

%%%%%%%%%%%%%%NORMALIZZAZIONE BLOCCHI%%%%%%%%%%%%%%%%%%

	[a1,s,t]=norm_blocco(a1);
if s==0
	a2=a2-s+1;
else
	a2=((a2-s)/t)+1;
end

	[b1,s,t]=norm_blocco(b1);
	
if s==0
	b2=b2-s+1;
else
	b2=((b2-s)/t)+1;
end
	[c1,s,t]=norm_blocco(c1);
if s==0
	c2=c2-s+1;
else
	c2=((c2-s)/t)+1;
end
	[d1,s,t]=norm_blocco(d1);
if s==0
	d2=d2-s+1;
else
	d2=((d2-s)/t)+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% media del quaternione q1
a1m=mean2(a1);
b1m=mean2(b1);
c1m=mean2(c1);
d1m=mean2(d1);

% media del quaternione q2
a2m=mean2(a2);
b2m=mean2(b2);
c2m=mean2(c2);
d2m=mean2(d2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_q1m=sqrt((a1m^2)+(b1m^2)+(c1m^2)+(d1m^2));
mod_q2m=sqrt((a2m^2)+(b2m^2)+(c2m^2)+(d2m^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod_q1=sqrt((a1.^2)+(b1.^2)+(c1.^2)+(d1.^2));
mod_q2=sqrt((a2.^2)+(b2.^2)+(c2.^2)+(d2.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% termine 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
termine2 = (mod_q1m*mod_q2m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% termine 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%
termine4 = ((mod_q1m^2)+(mod_q2m^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% termine 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
int1=(size1*size2)/((size1*size2)-1)*mean2(mod_q1.^2);
int2=(size1*size2)/((size1*size2)-1)*mean2(mod_q2.^2);
termine3=int1+int2-(size1*size2)/((size1*size2)-1)*(mod_q1m^2)-(size1*size2)/((size1*size2)-1)*(mod_q2m^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
sigma_z1 = sqrt(int1-(size1*size2)/((size1*size2)-1)*(mod_q1m^2));
sigma_z2 = sqrt(int2-(size1*size2)/((size1*size2)-1)*(mod_q2m^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% fattore mean_bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_bias = 2*termine2/termine4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(termine3==0)
	qs = mean_bias;
	qi = mean_bias;
	qj = mean_bias;
	qk = mean_bias;
else
	cbm = 2 / termine3;

%%%%%%%%%%%%%%%%%%% termine 1%%%%%%%%%%%%%%%%%%%
exs=(size1*size2)/((size1*size2)-1)*mean2((a1.*a2)+(b1.*b2)+(c1.*c2)+(d1.*d2));
exi=(size1*size2)/((size1*size2)-1)*mean2((-a1.*b2)+(b1.*a2)+(-c1.*d2)+(d1.*c2));
exj=(size1*size2)/((size1*size2)-1)*mean2((-a1.*c2)+(b1.*d2)+(c1.*a2)+(-d1.*b2));
exk=(size1*size2)/((size1*size2)-1)*mean2((-a1.*d2)+(-b1.*c2)+(c1.*b2)+(d1.*a2));

pros=(a1m*a2m)+(b1m*b2m)+(c1m*c2m)+(d1m*d2m);
proi=(-a1m*b2m)+(b1m*a2m)+(-c1m*d2m)+(d1m*c2m);
proj=(-a1m*c2m)+(b1m*d2m)+(c1m*a2m)+(-d1m*b2m);
prok=(-a1m*d2m)+(-b1m*c2m)+(c1m*b2m)+(d1m*a2m);

termine1s=exs-(size1*size2)/((size1*size2)-1)*pros;
termine1i=exi-(size1*size2)/((size1*size2)-1)*proi;
termine1j=exj-(size1*size2)/((size1*size2)-1)*proj;
termine1k=exk-(size1*size2)/((size1*size2)-1)*prok;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% simplified hypcpxCC ****
hypcpxCCs = termine1s;
hypcpxCCi = termine1i;
hypcpxCCj = termine1j;
hypcpxCCk = termine1k;

%display(strcat('CCB=',num2str(hypcpxCCs),' CCG=',num2str(hypcpxCCi),' CCR=',num2str(hypcpxCCj),' CCN=',num2str(hypcpxCCk),' Mean=',num2str(mean_bias),' CinContrast=',num2str(cbm)));
%%%%%%%%%%%%%%%%%universal image quality index%%%%%%%%%%%
qs = hypcpxCCs * mean_bias * cbm;
qi = hypcpxCCi * mean_bias * cbm;
qj = hypcpxCCj * mean_bias * cbm;
qk = hypcpxCCk * mean_bias * cbm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end