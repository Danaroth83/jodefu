function [Bands_to_sharpen,Bands_to_display]=load_Bands_to_sharpen(sensor,label,im_tag)

if nargin<=1, label='all'; end
if nargin<=2, im_tag=[]; end

if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'CHRIS'), sensor='CHR'; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'GE1'), sensor='GeoEye1'; end
if strncmpi(im_tag,'Beijing',7) && strcmpi(sensor,'WV3'), sensor='WV34bands'; end

switch sensor
    case 'HYP'
        %%% Band to display options
        Bands_to_display=[29,23,16]; % Wavelength 640, 579, 508 (suggested on user guide)
        % Bands_to_display=[28,23,16]; % Wavelength 630, 579, 508 (minimum SAM test with MS)
        % Bands_to_display=[32,23,16]; % Wavelength 671, 579, 508 (best match to MS-ALI)
        % Bands_to_display=[30,18,13]; % Wavelength 651, 529, 478 (middle of color spectrum)       
        % Bands_to_display=[27,22,16]; % Wavelength 620, 569, 508 (original toolbox)
        % Bands_to_display=[33,18,14]; % Wavelength 681, 529, 488 (lowest visual fidelty within overlap bands)
        % Bands_to_display=[32,20,14]; % Wavelength 671, 548, 488
        if strcmpi(label,'all')
            Bands_to_sharpen=1:242;
        elseif strcmpi(label,'unnoise')
            if strncmpi(im_tag,'Beijing',7)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220];
            elseif strncmpi(im_tag,'Sydney',6)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220]; %Unchecked
            elseif strncmpi(im_tag,'RdJ',3)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220]; %Unchecked
            elseif strncmpi(im_tag,'SaoPaulo',8)
                Bands_to_sharpen=[9:57,82:96,101:118,134:164,183,187:220];
            elseif strncmpi(im_tag,'SanFrancisco',12)
                Bands_to_sharpen=[11:57,82:96,101:119,134:164,183,188:189,193:215,217];
            elseif strncmpi(im_tag,'Sofia',5)
                Bands_to_sharpen=[11:57,82:96,101:119,134:164,183:184,187:220];
            elseif strncmpi(im_tag,'Sudbury',7)
                Bands_to_sharpen=[11:57,82:97,99:119,133:164,188:189,193:218];
            else
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220];
            end
        elseif strcmpi(label,'VNIR')
            Bands_to_sharpen=load_Bands_to_sharpen(sensor,'unnoise',im_tag);
            Bands_to_sharpen=intersect(Bands_to_sharpen,9:57);
        elseif strcmpi(label,'16')
            Bands_to_sharpen=[14:16,18:24,28:33];
        elseif strcmpi(label,'MSovlp')
            [~,sensor_ovl]=load_default_sensor(place);
            Bands_to_sharpen=load_Band_overlap(sensor,sensor_ovl,'HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'16',im_tag));
        elseif strcmpi(label,'MSovlp2')
            Bands_to_sharpen=load_Band_overlap(sensor,'ALI','HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'16',im_tag));
        elseif strcmpi(label,'MSovlpVNIR')
            [~,sensor_ovl]=load_default_sensor(place);
            Bands_to_sharpen=load_Band_overlap(sensor,sensor_ovl,'HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'VNIR',im_tag));
        elseif strcmpi(label,'MSovlp2VNIR')
            Bands_to_sharpen=load_Band_overlap(sensor,'ALI','HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'VNIR',im_tag));
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[16,23,29];
        else
            Bands_to_sharpen=14:33;
        end
    case 'ALI'
        Bands_to_display=2:4;
        if any(strcmpi(label,{'all','unnoise'}))
            Bands_to_sharpen=1:9;
        elseif strcmpi(label,'VNIR')
            Bands_to_sharpen=1:6;
        elseif strcmp(label,'16')
            Bands_to_sharpen=1:4;
        elseif strcmp(label,'RGB')
            Bands_to_sharpen=2:4;
        else
            Bands_to_sharpen=1:9;
        end
    case {'WV2','WV3'}
        if any(strcmpi(label,{'RGB','3'}))
            Bands_to_sharpen=[2,3,5];
        elseif strcmpi(label,'4')
            Bands_to_sharpen=[2,3,5,7];
        else
            Bands_to_sharpen=1:8;
        end
        Bands_to_display=[5,3,2];
    case {'WV24bands','WV34bands'}
        if any(strcmpi(label,{'RGB','3'}))
            Bands_to_sharpen=[1,2,3];
        else
            Bands_to_sharpen=1:4;
        end
        Bands_to_display=[3,2,1];
    case 'ROSIS'
        if strcmpi(label,'32')
            Bands_to_sharpen=12:2:74; %32
            Bands_to_display=[60,30,16]; % Wavelength 491, 547, 667
        elseif strcmpi(label,'8')
            Bands_to_sharpen=10:10:80; %8
            Bands_to_display=[60,30,20];
        elseif strcmpi(label,'16')
            Bands_to_sharpen=10:5:85; %16
            Bands_to_display=[60,30,15];
        elseif strcmpi(label,'last')
            Bands_to_sharpen=51:1:66; %16
            Bands_to_display=[60,56,51];
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[16,30,60];
            Bands_to_display=[60,30,16]; % Wavelength 491, 547, 667
        else
            Bands_to_sharpen=1:103;
            Bands_to_display=[60,30,16];
        end
    case {'AVIRIS','AVIRIS220','AVIRIS176'}
        Bands_to_display=[21,9,3]; % Bands ?? watch out: to choose
        if strcmpi(label,'all')
            Bands_to_sharpen=1:176; % All bands
        elseif strcmpi(label,'4')
            Bands_to_sharpen=[3,9,15,21];
        elseif strcmpi(label,'ovlp')
            Bands_to_sharpen=1:41;
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[3,9,21];
        else
            Bands_to_sharpen=1:2:32;
        end
    case 'CHR'
        if strcmpi(label,'16')
            Bands_to_sharpen=3:17; 
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[7,9,11];
        else
            Bands_to_sharpen=1:18;
        end
        Bands_to_display=[11,9,7];
    case 'DE2'
        Bands_to_sharpen=1:4;
        Bands_to_display=[4,3,2];
    case {'CAVE','U260'}
        if strcmpi(label,'16')
            Bands_to_sharpen=1:2:31;
        elseif strcmpi(label,'4')
            Bands_to_sharpen=7:8:31;
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[7,16,22];
        else
            Bands_to_sharpen=1:31;
        end
        % Bands_to_display=[22,16,6]; %wavelengths: [610,550,450] nm
        % Bands_to_display=[31,16,10]; %wavelengths: [700,550,490] nm
        % Bands_to_display=[29,15,9]; %wavelengths: [680,540,480] nm
        % Bands_to_display=[28,11,5]; %wavelengths: [670,500,440] nm
        Bands_to_display=[22,16,7]; %wavelengths: [610,550,460] nm
    case {'CZ','Nuance','TT','Varispec'}
        if strcmpi(label,'16')
            Bands_to_sharpen=1:2:31;
        elseif strcmpi(label,'4')
            Bands_to_sharpen=3:8:27;
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[5,14,20];
        else
            Bands_to_sharpen=1:31;
        end
        % Bands_to_display=[20,14,4]; %wavelengths: [610,550,450] nm
        % Bands_to_display=[29,12,6]; %wavelengths: [700,530,470] nm
        % Bands_to_display=[27,13,7]; %wavelengths: [680,540,480] nm
        % Bands_to_display=[26,10,3]; %wavelengths: [670,510,440] nm
        Bands_to_display=[20,14,5]; %wavelengths: [610,550,460] nm
    case {'FNA','Pulnix'}
        if strcmpi(label,'16')
           Bands_to_sharpen=2:2:32;
        elseif strcmpi(label,'4')
            Bands_to_sharpen=4:8:28;
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[6,15,21];
        else
            if strncmpi(im_tag,'FNAc',4)
                if strncmpi(im_tag,'FNAh',4)
                    Bands_to_sharpen=1:32; 
                else
                    Bands_to_sharpen=1:33;
                end                
            elseif strncmpi(im_tag,'FNAc',4)
                Bands_to_sharpen=1:33;
            else
                Bands_to_sharpen=2:32;
            end
        end
        Bands_to_display=[21,15,6]; %wavelengths: [610,550,460] nm
    otherwise
        if any(strcmpi(label,{'RGB','3'}))
            Bands_to_sharpen=1:3;
        else
            Bands_to_sharpen=1:4;
        end
        if any(strcmpi(sensor,{'WV3','WV2','QB','IKONOS','GeoEye1','Pleiades'}))
            Bands_to_display=3:-1:1;
        else
            Bands_to_display=1:3;
        end
end
        
end
            