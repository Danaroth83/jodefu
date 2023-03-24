function [mse_index,mse_map] = MSE( I_1,I_2 )

nelem=numel(I_1);
if nelem~=numel(I_2), error('Dimensions of compared images must be equal'); end

mse_map=(I_1-I_2).^2;
mse_index=sum(mse_map(:))/nelem;


end

