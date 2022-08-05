function I_out = interp23tapGeneral_mod3(I_in,ratio,L,window,flag_edge)
% keyboard
if nargin<=2 || isempty(L), L = 44; end % tap
if nargin<=3 || isempty(window), window='Hamming'; end
if nargin<=4 || isempty(flag_edge), flag_edge='circular'; end % Also available 'symmetric'
if numel(ratio)==1, ratio=[ratio,ratio]; end
if numel(L)==1, L=[L,L]; end

L=L+(rem(ratio,2)==rem(L,2));

switch window
    case 'Rectangular'
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),rectwin(L(1)+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),rectwin(L(2)+1),'noscale');
    case 'Blackman'
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),blackman(L(1)+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),blackman(L(2)+1),'noscale');
    otherwise
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),'noscale');
end

I_out=zeros(size(I_in,1)*ratio(1),size(I_in,2)*ratio(2),size(I_in,3));
I_in=repelem(I_in,2-rem(ratio(1),2),2-rem(ratio(2),2));
idx=sort(unique(cat(2,ceil(ratio(1)/2):ratio(1):size(I_out,1),floor(ratio(1)/2)+1:ratio(1):size(I_out,1))));
idy=sort(unique(cat(2,ceil(ratio(2)/2):ratio(2):size(I_out,2),floor(ratio(2)/2)+1:ratio(2):size(I_out,2))));
I_out(idx,idy,:)=I_in;

for ii = 1 : size(I_in,3)
    t = I_out(:,:,ii);
    if rem(ratio(1),2)==0
        t1=t(1:2:end,:,:);
        t1=imfilter(t1',BaseCoeff_x(1:2:end),flag_edge);
        t2=t(2:2:end,:,:);
        t2=imfilter(t2',BaseCoeff_x(2:2:end),flag_edge);
        t(1:2:end,:,:)=t1';
        t(2:2:end,:,:)=t2';
    else
        t = (imfilter(t',BaseCoeff_x,flag_edge))';
    end
    if rem(ratio(2),2)==0
        t1=t(:,1:2:end,:);
        t1=imfilter(t1,BaseCoeff_y(1:2:end),flag_edge);
        t2=t(:,2:2:end,:);
        t2=imfilter(t2,BaseCoeff_y(2:2:end),flag_edge);
        I_out(:,1:2:end,ii)=t1;
        I_out(:,2:2:end,ii)=t2;
    else
        I_out(:,:,ii) = imfilter(t,BaseCoeff_y,flag_edge);
    end
end

end