function I_PAN_LR = LPF_filter(I_PAN,Resize_fact,LR_method,downsample,GNyqPan)

if nargin<=1 || isempty(LR_method)
    LR_method='resize';
    % LR_method='MTF';
    % LR_method='filter';
end
if nargin<=2 || isempty(downsample)
    downsample=true;
end
switch LR_method
    case {1,'imresize'}
        LR_method='resize';
    case {0,'MTFfiltered'}
        LR_method='MTF';
end

%% Low Resolution PAN construction
switch LR_method
    case 'resize'
        if downsample
            I_PAN_LR=imresize(I_PAN,1/Resize_fact);
        else
            I_PAN_LR=imresize(imresize(I_PAN,1/Resize_fact),Resize_fact);
        end
    case 'filter'
        fcut=1/Resize_fact;
        N = 21;
        [f1,f2] = freqspace(21,'meshgrid');
        Hd = ones(N);
        Hd((abs(f1)>fcut)|(abs(f2)>fcut)) = 0;
        h = fwind1(Hd,hamming(N));
        I_PAN_LR = imfilter(I_PAN,real(h),'replicate');
        if downsample
            I_PAN_LR = imresize(I_PAN_LR,1/Resize_fact,'nearest');
        end
    case 'MTF'
        % switch sensor
        %     case 'QB'
        %         GNyqPan = 0.15;
        %     case 'IKONOS'
        %         GNyqPan = 0.17;
        %     case 'GeoEye1'
        %         GNyqPan = 0.16;
        %     case 'WV2'
        %         GNyqPan = 0.11;
        % end
        if nargin<=3
            GNyqPan=0.15;
        end
        
        N = 41;
        fcut = 1 / Resize_fact;
        alpha = sqrt(((N-1)*(fcut/2))^2/(-2*log(GNyqPan)));
        H = fspecial('gaussian', N, alpha);
        Hd = H./max(H(:));
        h = fwind1(Hd,kaiser(N));
        I_PAN_LR = imfilter(I_PAN,real(h),'replicate');
        if downsample
            I_PAN_LR = imresize(I_PAN_LR,1/Resize_fact,'nearest');
        end
end