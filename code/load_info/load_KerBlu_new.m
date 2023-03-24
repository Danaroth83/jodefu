function  [MTF,KerBlu,start_pos] = load_KerBlue_new( sensor, type, bands,ratio )
%LOAD_MTF Loads the Nyquist gain of the MS sensor MTF

if nargin<=0, sensor='none'; end
if nargin<=1, ratio=4; end
if nargin<=2, type='PAN'; end
if nargin<=4, im_tag=[]; end

if strcmpi(type,'PAN')
    switch sensor
        case 'IKONOS'
            MTF=0.17;
        case 'GeoEye1'
            MTF=0.16;
        case 'QB'
            MTF=0.15;
        case 'WV2'
            MTF=0.11;
        case {'WV3','WV34bands'}
            MTF=0.14;
        case 'ALI'
            MTF=0.13;
        otherwise
            MTF=0.15;
    end
else
    switch sensor
        case 'IKONOS'
            MTF=[0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
        case 'GeoEye1'
            MTF=[0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
        case 'QB'
            MTF=[0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
        case 'WV2'
            MTF=[0.35 .* ones(1,7), 0.27];
        case 'WV2_test'
            MTF = 0.15 .* ones(1,8);
        case 'WV3'   
            MTF=[0.32,0.36,0.36,0.35,0.36,0.36,0.33,0.32]; % Band Order: Coastal,B,G,Y,R,R-edge,NIR1,NIR2
        case 'WV34bands'
            MTF=[0.36 .* ones(1,3), 0.33]; % 4 bands bundle uses bands 2,3,5,7 (B,G,R,NIR1)
        case 'ALI'
            MTF=[0.29,0.30,0.28,0.29,0.28,0.29,0.25,0.25,0.25];
        case 'HYP'
            % MTF_MS=[0.27*ones(1,70), 0.29*ones(1,172)];
            %VNIR
            MTF(1:21)=0.27;
            MTF(22:41)=0.28;
            MTF(42:49)=0.26;
            MTF(50:70)=0.26;
            %SWIR
            MTF(71:100)=0.30;
            MTF(101:130)=0.30;
            MTF(131:177)=0.27;
            MTF(177:242)=0.27;
            % Error_Measure=0.03;
            % MTF_MS=MTF_MS+Error_Measure;
        case 'AVIRIS'
            MTF=0.45 .* ones(1,224);
        case 'AVIRIS_176'
            MTF=0.45 .* ones(1,176);
        otherwise
            if nargin==3
                MTF = 0.29 .* ones(1,length(bands));
                bands=1:length(bands);
                flag_failedloadingMTF=1;
            else
                error('Error: Please specify either bands or known sensors');
            end
    end

    if flag_failedloadingMTF==0
        if isempty(bands)
            MTF=NaN;
        elseif max(bands)<=length(MTF) && min(bands)>=1
            MTF=MTF(bands);
        else
            error('Specified bands are outside of known range');
        end
    end
end

if nargout >=2
    

end

