function h = genMTF(Resize_fact, sensor, varargin)
% Generate a bank of filters shaped on the MTF of the sensor. Each filter
% corresponds to a band acquired by the sensor. The 
% Examples
%   h = genMTF(4, 'QB');
%   h = genMTF(4, 'none', 'WV2');


switch sensor
    case 'QB'
        GNyq = [0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
    case 'IKONOS'
        GNyq = [0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
    case 'GeoEye1'
        GNyq = [0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
    case 'WV2'
        GNyq = [0.35 .* ones(1,7), 0.27];
    case {'WV3','WV3_4bands'}
        GNyq = [0.325 0.355 0.360 0.350 0.365 0.360 0.335 0.315];
        if strcmp(sensor,'WV3_4bands')
            GNyq = GNyq([2 3 5 7]);
        end
    case {'HYP','HYP_14_33','HYP_16_31'}
        flag_resize_new = 2; % MTF usage
        %VNIR
        MTF_MS(1:21)=0.27;
        MTF_MS(22:41)=0.28;
        MTF_MS(42:49)=0.26;
        MTF_MS(50:70)=0.26;
        %SWIR
        MTF_MS(71:100)=0.30;
        MTF_MS(101:130)=0.30;
        MTF_MS(131:177)=0.27;
        MTF_MS(177:242)=0.27;
        if strcmp(sensor,'HYP_14_33')
            GNyq = MTF_MS(14:33);
        elseif strcmp(sensor,'HYP_16_31')
            GNyq = MTF_MS(16:31);
        else
            GNyq = MTF_MS;
        end
    case 'none'
        GNyq = 0.29 .* ones(1,size(I_MS,3));
end


%%% MTF
N = 41;
nBands = length(GNyq);
h = zeros(N, N, nBands);
fcut = 1/Resize_fact;
   
for ii = 1 : nBands
    % Shape finestra gaussiana (Hp: in Fast and Efficient Panchromatic
    % Sharpening, TGARS 2010)
    alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h(:,:,ii) = fwind1(Hd,kaiser(N));
end