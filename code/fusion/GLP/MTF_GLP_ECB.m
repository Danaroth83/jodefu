%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           MTF_GLP_ECB fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP) with ECB algorithm. 
% 
% Interface:
%           I_Fus_MTF_GLP_ECB = MTF_GLP_ECB(I_MS,I_PAN,ratio,blockSize,c,GNyq)
%
% Inputs:
%           I_MS:               MS image upsampled at PAN scale;
%           I_PAN:              PAN image;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value;
%           blockSize:          Block size (e.g. 9 x 9);
%           c:                  Parameter between 2 and 3 (e.g. 2.5);
%           GNyq:               MTF gain of MS sensor at Nyquist frequency
%
% Outputs:
%           I_Fus_MTF_GLP_ECB:  MTF_GLP_ECB pansharpened image.
% 
% References:
%           [Aiazzi06]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, �MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,�
%                       Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591�596, May 2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_Fus_MTF_GLP_ECB = MTF_GLP_ECB(I_MS,I_PAN,ratio,blockSize,c,GNyq)

%%% Check BlockSize
if rem(blockSize,2) == 0
    %%% Approximate to odd BlockSize
    blockSize = blockSize + 1; 
end


imageHR = double(I_PAN);
I_MS = double(I_MS);

imageHR = repmat(imageHR,[1 1 size(I_MS,3)]);

% switch sensor
%     case 'QB' 
%         GNyq = [0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
%     case 'IKONOS'
%         GNyq = [0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
%     case 'GeoEye1'
%         GNyq = [0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
%     case 'WV2'
%         GNyq = [0.35 .* ones(1,7), 0.27];
%     case 'none'
%         if strcmp(tag,'WV2')
%             GNyq = 0.15 .* ones(1,8);
%         else
%             GNyq = 0.29 .* ones(1,size(I_MS,3));
%         end
% end


%%% MTF
N = 41;
PAN_LP = zeros(size(I_MS));
nBands = size(I_MS,3);
fcut = 1/ratio;
   
for ii = 1 : nBands
    alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    PAN_LP(:,:,ii) = imfilter(imageHR(:,:,ii),real(h),'replicate');
end

PAN_LP = double(PAN_LP);
PAN_LP = imresize(PAN_LP,1/ratio,'nearest');
if (2^round(log2(ratio)) ~= ratio)
    imageHRLP = imresize(PAN_LP,ratio);
else
    imageHRLP = interp23tap(PAN_LP,ratio);
end

DetailsHRPan = imageHR - imageHRLP;

I_Fus_MTF_GLP_ECB = zeros(size(I_MS));
for b = 1 : size(I_MS,3)

     Global_Correlation = corr2(I_MS(:,:,b),imageHRLP(:,:,b));

     for y=1 : size(I_MS,1)
        for x=1 : size(I_MS,2)
            
            startx = x - floor(blockSize/2);
            starty = y - floor(blockSize/2);
            endy = y + floor(blockSize/2);
            endx = x + floor(blockSize/2);
            
            if endy > size(I_MS,1)
                endy = size(I_MS,1);
            end
            
            if endx > size(I_MS,2)
                endx = size(I_MS,2);
            end
            
            if starty < 1
                starty = 1;
            end
            
            if startx < 1
                startx = 1;
            end
            
            BlockMS = I_MS((starty:endy),(startx:endx),b);
            BlockPan = imageHRLP((starty:endy),(startx:endx),b);
            StdMS=std2(BlockMS);
            StdPan=std2(BlockPan);
            Correlation=corr2(BlockMS,BlockPan);

            LG = min((Correlation./Global_Correlation) .* (StdMS/StdPan),c);
            
            BlockDetails = DetailsHRPan(y,x,b);
            Details2Add=LG*BlockDetails;
            FusedBlock = I_MS(y,x,b) + Details2Add;
            I_Fus_MTF_GLP_ECB(y,x,b)=FusedBlock;
        end
    end
end

end