function h = genPANMTF(Resize_fact, sensor, varargin)
% Generate a bank of filters shaped on the MTF of the sensor. Each filter
% corresponds to a band acquired by the sensor. The
% Examples
%   h = genMTF(4, 'QB');
%   h = genMTF(4, 'none', 'WV2');


switch sensor
    case 'QB'
        GNyq = 0.15;
    case 'IKONOS'
        GNyq = 0.17;
    case 'GeoEye1'
        GNyq = 0.16;
    case 'WV2'
        GNyq = 0.11;
    case {'WV3','WV3_4bands'}
        GNyq = 0.14;
    case {'HYP','HYP_14_33','HYP_16_31'}
        GNyq = 0.11;
    case 'none'
        GNyq = 0.15;
end

%%% MTF
N = 41;
fcut = 1/Resize_fact;

alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq)));
H = fspecial('gaussian', N, alpha);
Hd = H./max(H(:));
h = fwind1(Hd,kaiser(N));
