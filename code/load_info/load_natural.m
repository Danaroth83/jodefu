%% Loading natural image
%
% Author:
% Daniele Picone
% 
% Description:
% Routine for loading data and informations relative to natural datasets
%
% Usage:
% I_load=load_natural(im_tag)
%
% Input:
% im_tag: String describing the input dataset, in the form
%         [Dataset_label][Image_number]_[cut_label]_[Band_label]
%         eg: 'CAVE4_cut2_16'
%         Available Dataset_labels: 'CAVE','FNAc','FNAh','FNAl','CZi','CZo'
%         Available cut_label: 'cutall','cut1','cut2','cut3' (default: 'cutall')
%         Available Band_label: '4','16','RGB','all' (default: 'all')
%
% Output:
% I_load: Struct containing data and information of the loaded dataset,
%         whose fields are:
%      - data: The image itself
%      - sensor: String describing the sensor
%      - im_tag: the string describing the dataset
%      - Bands_selected: List of bands loaded from the full dataset
%      - Bands_to_display: List of 3 representative bands to show
%      - DynamicRange: Maximum possible radiometric value
%      - IntBits: Bits used to represent the dataset
%      - wavelength: Central wavelength of each loaded channel
%      - bandwidth: Bandwidth of each loaded channel
%
%
% References:
% [1] F. Yasuma, T. Mitsunaga, D. Iso, and S.K. Nayar. "Generalized 
%     Assorted Pixel Camera: Post-Capture Control of Resolution, Dynamic 
%     Range and Spectrum," Technical Report, Department of Computer 
%     Science, Columbia University CUCS-061-08, Nov. 2008.
% [2] Chakrabarti A. and Zickler T. "Statistics of Real-World Hyperspectral
%     Images," in Proceedings of the IEEE Conference on CVPR, 2011
% [3] Nascimento, S.M.C., Ferreira, F.P., & Foster, D.H. (2002). 
%     "Statistics of spatial cone-excitation ratios in natural scenes." 
%     Journal of the Optical Society of America A, 19, 1484-1490.


function I_load=load_natural(im_tag)

current_folder=fileparts(mfilename('fullpath'));
data_folder=fullfile(current_folder,'..','..','data','Natural');

%% Parse input
[place,cut,Bands,sensor]=parse_im_tag(im_tag);
type='HS';
[vcut,hcut,edge_cut,~,Qblocks_size]=load_cut(place,cut,sensor);
[Bands_to_sharpen,Bands_to_display]=load_Bands_to_sharpen(sensor,Bands,place);
[~,Bands_to_display]=min(abs(repmat(Bands_to_sharpen',[1,3])-repmat(Bands_to_display,[length(Bands_to_sharpen),1])),[],1);
[GSD_GT,L]=load_resolution(sensor, im_tag, type);
[wl_central,wl_bandwidth]=load_wavelength(sensor,'HS',Bands_to_sharpen);

%% Get image index number
im_num_new=1; tail=0;
while ~isnan(im_num_new)
    im_num=im_num_new;
    im_num_new=str2double(place(end-tail:end)); 
    tail=tail+1; 
end


if strncmpi(place,'CZ',2)
    if strcmpi(place(3),'i')
        file_folder=fullfile('CZ','CZ_hsdbi');
        filelist={...
            'img3.mat' ,'img4.mat' ,'img5.mat' ,'img6.mat' ,'imga4.mat',...
            'imga8.mat','imgc3.mat','imgc6.mat','imgd0.mat','imgd1.mat',...
            'imgd5.mat','imgd6.mat','imgg0.mat','imgg1.mat','imgg2.mat',...
            'imgg3.mat','imgg4.mat','imgg5.mat','imgg6.mat','imgg7.mat',...
            'imgg8.mat','imgg9.mat','imgh4.mat','imgh5.mat','imgh6.mat',...
            'imgh6.mat'};
    else
        file_folder=fullfile('CZ','CZ_hsdb');
        filelist={...
            'img1.mat' ,'img2.mat' ,'imga1.mat','imga2.mat','imga5.mat',...
            'imga6.mat','imga7.mat','imgb0.mat','imgb1.mat','imgb2.mat',...
            'imgb3.mat','imgb4.mat','imgb5.mat','imgb6.mat','imgb7.mat',...
            'imgb8.mat','imgb9.mat','imgc1.mat','imgc2.mat','imgc4.mat',...
            'imgc5.mat','imgc6.mat','imgc7.mat','imgc8.mat','imgc9.mat',...
            'imgd2.mat','imgd3.mat','imgd4.mat','imgd7.mat','imgd8.mat',...
            'imgd9.mat','imge0.mat','imge1.mat','imge2.mat','imge3.mat',...
            'imge4.mat','imge5.mat','imge6.mat','imge7.mat','imgf1.mat',...
            'imgf2.mat','imgf3.mat','imgf4.mat','imgf5.mat','imgf6.mat',...
            'imgf7.mat','imgf8.mat','imgh0.mat','imgh1.mat','imgh2.mat',...
            'imgh3.mat'};
    end
    load(fullfile(data_folder,file_folder,filelist{im_num}),'ref');
    I_in = double(ref(vcut,hcut,Bands_to_sharpen));
    fid  = fopen(fullfile(data_folder,file_folder,'calib.txt'));
    I_load.signature = fscanf(fid,'%f');
    fclose(fid);
    I_in=I_in ./ shiftdim(I_load.signature(Bands_to_sharpen),-2);
elseif strncmpi(place,'CAVE',4)
    file_folder=fullfile('CAVE');
    filelist={...
        'balloons','beads','cd','chart_and_stuffed_toy','clay',...
        'cloth','egyptian_statue','face','fake_and_real_beer',...
        'fake_and_real_food','fake_and_real_lemon_slices',...
        'fake_and_real_lemons','fake_and_real_peppers',...
        'fake_and_real_strawberries.','fake_and_real_sushi',...
        'fake_and_real_tomatoes','feathers','flowers',...
        'glass_tiles','hairs','jelly_beans','oil_painting.',...
        'paints','photo_and_face','pompoms.','real_and_fake_apples',...
        'real_and_fake_peppers','sponges','stuffed_toys.',...
        'superballs','thread_spools','watercolors'};
    load(fullfile(data_folder,file_folder,[filelist{im_num},'_ms.mat']),'I_MS');
    I_in = double(I_MS(vcut,hcut,Bands_to_sharpen));
elseif strncmpi(place,'FNA',3)
    if strcmpi(place(4),'h')
        file_folder=fullfile('FNA','HR',sprintf('scene%d',im_num));
        filelist={...
            'ref_crown3bb_reg1_lax.mat','ref_ruivaes1bb_reg1_lax.mat',...
            'ref_mosteiro4bb_reg1_lax.mat','ref_cyflower1bb_reg1.mat',...
            'ref_cbrufefields1bb_reg1_lax.mat','ref_braga1bb_reg1.mat',...
            'ref_ribeira1bbb_reg1.mat','ref_farme1bbbb_reg1.mat'};
        load(fullfile(data_folder,file_folder,filelist{im_num}),'reflectances');
        I_in = double(imcrop_custom(reflectances(:,:,Bands_to_sharpen),vcut,hcut));
    elseif strcmpi(place(4),'c')
        file_folder=fullfile('FNA','Closeup');
        filelist={...
            'Bom_Jesus_Bush','Bom_Jesus_Marigold','Bom_Jesus_Red_flower',...
            'Bom_Jesus_Ruin','Braga_Grafitti','Gualtar_Columns',...
            'Gualtar_Orange_Trees','Gualtar_Steps','Gualtar_Villa',...
            'Lillies_Closeup','Lilly_Closeup','Ruivaes_Fern',...
            'Ruivaes_Ferns','Ruivaes_Ruin','Sameiro_Bark',...
            'Sameiro_Forest','Sameiro_Glade','Sameiro_Leaves',...
            'Sameiro_Trees','Sete_Fontes_Rock','Souto_Farm_Barn',...
            'Souto_Farm_Barn','Souto_Roog_Tiles','Souto_Wood_Pile',...
            'Tenoes_Wall','Tenoes_Wall_Closeup','Tibaes_Corridor',...
            'Tibaes_Garden','Tibaes_Garden_Entrance','Yellow_Rose'};
        load(fullfile(data_folder,file_folder,filelist{im_num},[filelist{im_num},'.mat']),'hsi');
        I_in = double(imcrop_custom(hsi(:,:,Bands_to_sharpen),vcut,hcut));
    else
        file_folder=fullfile('FNA','LR',sprintf('scene%d',im_num));
        load(fullfile(data_folder,file_folder,sprintf('scene%d.mat',im_num)),'reflectances');
        I_in = double(imcrop_custom(reflectances(:,:,Bands_to_sharpen-1),vcut,hcut)); % Note: The -1 is an hack to be consistent with the definition of Bands_to_sharpen
    end
elseif strncmpi(place,'TT',2)
    file_folder=fullfile('TT');
    filelist={...
            'Butterfly','Butterfly2','Butterfly3','Butterfly4',...
            'Butterfly5','Butterfly6','Butterfly7','Butterfly8',...
            'CD','Character','Chart24','ChartCZP','ChartDC',...
            'ChartRes','ChartSG','Cloth','Cloth2','Cloth3','Cloth4',...
            'Cloth5','Cloth6','Color','Colorchart','Doll','Fan',...
            'Fan2','Fan3','Flower','Flower2','Flower3','Party',...
            'Tape','Tape2','Tshirts','Tshirts2'};
    load(fullfile(data_folder,file_folder,[filelist{im_num},'.mat']),'img');
    I_in = double(imcrop_custom(img(:,:,Bands_to_sharpen),vcut,hcut));
elseif strncmpi(place,'Moffett',7)
    file_folder=fullfile('..','HS');
    load(fullfile(data_folder,file_folder,'Moffett_HS_RR.mat'),'I_MS');
    I_in = double(imcrop_custom(I_MS(:,:,Bands_to_sharpen),vcut,hcut));
elseif strncmpi(place,'Pavia',5)
    file_folder=fullfile('..','HS');
    load(fullfile(data_folder,file_folder,'Pavia_RR.mat'),'I_MS');
    I_in = double(imcrop_custom(I_MS(:,:,Bands_to_sharpen),vcut,hcut));
end

I_load.data             = I_in;
I_load.sensor           = sensor;
I_load.type             = type;
I_load.im_tag           = im_tag;
I_load.Bands_selected   = Bands_to_sharpen;
I_load.Bands_to_display = Bands_to_display;
I_load.DynamicRange     = max(2^L-1,1);
I_load.IntBits          = L;
I_load.edge_cut         = edge_cut;
I_load.ratio            = 1;
I_load.Qblocks_size     = Qblocks_size;
I_load.GSD              = GSD_GT;
I_load.label            = 'GT';
I_load.wavelength       = wl_central;
I_load.bandwidth        = wl_bandwidth;
I_load.size             = [size(I_load.data,1),size(I_load.data,2),size(I_load.data,3)];
% I_load.signature        = signature;

% I_req{ii}.Band_overlap_PAN=Band_overlap_PAN;
% I_req{ii}.GNyq=GNyq;
% I_req{ii}.comp_factor=comp_factor;
% I_req{ii}.comp_method=comp_method;
% I_req{ii}.spectralweights=spectralweights;
