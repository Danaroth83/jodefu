function I_out = interp23tap_mod4(I_in,ratio,order,flag_edge)
% keyboard
if numel(ratio)==1, ratio=[ratio,ratio]; end
if nargin<=2 || isempty(order), order = 11; end % Lagrange polynomial interpolation order
if nargin<=3, flag_edge='circular'; end         % Also available 'symmetric','replicate'

if rem(ratio(1),2)==0,
    BaseCoeff_x = intfilt(ratio(1)*2,order,'Lagrange');
    BaseCoeff_x = BaseCoeff_x(rem(order,2)+1:2:end);
else
    BaseCoeff_x = intfilt(ratio(1),order,'Lagrange');
end
if rem(ratio(2),2)==0,
    BaseCoeff_y = intfilt(ratio(2)*2,order,'Lagrange');
    BaseCoeff_y = BaseCoeff_y(rem(order,2)+1:2:end);
else
    BaseCoeff_y = intfilt(ratio(2),order,'Lagrange');
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