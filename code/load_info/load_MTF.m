function  MTF_MS = load_MTF( sensor_MS, im_tag, bands, type )
%LOAD_MTF Loads the Nyquist gain of the MS sensor MTF

if nargin<=1, im_tag=[]; end
if nargin<=3, type='MS'; end
if strcmpi(type,'PAN')
    MTF_MS=load_MTF_PAN(sensor_MS);
else
    flag_failedloadingMTF=0;

    switch sensor_MS
        case 'IKONOS'
            MTF_MS=[0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
        case 'GeoEye1'
            MTF_MS=[0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
        case 'QB'
            MTF_MS=[0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
        case 'WV2'
            MTF_MS=[0.35 .* ones(1,7), 0.27];
        case 'WV3'
            if strncmpi(im_tag,'Beijing',7)
                MTF_MS=[0.36 .* ones(1,3), 0.33]; % 4 bands bundle uses bands 2,3,5,7 (B,G,R,NIR1)
            else    
                MTF_MS=[0.32,0.36,0.36,0.35,0.36,0.36,0.33,0.32]; % Band Order: Coastal,B,G,Y,R,R-edge,NIR1,NIR2
            end
        case 'WV34bands'
            MTF_MS=[0.36 .* ones(1,3), 0.33]; % 4 bands bundle uses bands 2,3,5,7 (B,G,R,NIR1)
        case 'ALI'
            MTF_MS=[0.29,0.30,0.28,0.29,0.28,0.29,0.25,0.25,0.25];
        case 'HYP'
            % MTF_MS=[0.27*ones(1,70), 0.29*ones(1,172)];
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
            % Error_Measure=0.03;
            % MTF_MS=MTF_MS+Error_Measure;
        case 'AVIRIS'
            if strncmpi(im_tag,'Moffett',7)
                MTF_MS=0.45 .* ones(1,176);
            else
                MTF_MS=0.45 .* ones(1,224);
            end
        otherwise
            if strncmpi(im_tag,'WV2',3)
                MTF_MS = 0.15 .* ones(1,8);
            elseif nargin==3
                MTF_MS = 0.29 .* ones(1,length(bands));
                bands=1:length(bands);
                flag_failedloadingMTF=1;
            else
                error('Error: Please specify either bands or known sensors');
            end
    end

    if flag_failedloadingMTF==0
        if isempty(bands)
            MTF_MS=NaN;
        elseif max(bands)<=length(MTF_MS) && min(bands)>=1
            MTF_MS=MTF_MS(bands);
        else
            error('Specified bands are outside of known range');
        end
    end
end

end

