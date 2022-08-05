function [I_out, I_GT, MR] = wrapper_software_compression(varargin)

%% Support folders' paths
current_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(current_folder,'..','..','data','output');
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'load_info'),...
        fullfile(project_folder,'demosaic'),...
        fullfile(project_folder,'fusion'),...
        fullfile(project_folder,'interpolation'),...
        fullfile(project_folder,'mosaic'),...
        fullfile(project_folder,'validation'),...
        fullfile(project_folder,'visualization'));

im_tag = 'Washington_cut256_RGB';
compression_ratio = 0.25;
qindex_list = {'SSIM','PSNR','ERGAS','SAM','sCC','UIQI','Q2n'};
flag_vis = 1;
methods_list = {'BIN', 'JPEG'};
output_fol = "software_compression";


for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1}; 
    if any(strcmpi(pname,{'im_tag','image','img','im'}))                        % Image tag
        im_tag = pval;
    elseif any(strcmpi(pname,{'compression','compression_ratio'}))              % Target compression ratio
        compression_ratio = pval;
    elseif any(strcmpi(pname,{'qindex_list','qindex','index','index_list'}))    % List of quality indices for validation
        qindex_list = pval;
    elseif any(strcmpi(pname,{'figure','visualization','vis'}))                 % Flag for visualizing images
        flag_vis = pval;
    elseif any(strcmpi(pname,{'method','method_list'}))                         % Software compression method
        methods_list = pval;
    elseif any(strcmpi(pname,{'output','folder','output_folder','output_fol'})) % Output folder
        output_fol=pval;
    end
end
output_folder = fullfile(output_folder, output_fol);
I_load=load_dataset_pansharpening(im_tag,'request',{'GT'});
I_GT = I_load{1};

reference = I_GT.data;
I_out_data = zeros(size(reference,1), size(reference,2), size(reference,3), numel(methods_list));
for kk = 1:numel(methods_list)
    method_current = methods_list{kk};
    dynamic_range = I_GT.DynamicRange;
    reference = I_GT.data;
    if any(strcmpi(method_current,{'bin','binning'}))
        dynamic_range_new = ceil(2^(log2(dynamic_range+1) * compression_ratio))-1;
        I_out_data(:,:,:,kk) = round(reference / dynamic_range * dynamic_range_new) * dynamic_range / dynamic_range_new;
    elseif any(strcmpi(method_current,{'jpeg'}))
        dynamic_range_new = 255;
        I_temp = round(reference / dynamic_range * dynamic_range_new);
        new_compression_ratio = 1 / compression_ratio ./ log2(dynamic_range+1) * 8;
        file_temp = 'temp.jp2';
        for mm = 1:size(I_temp,3)
            imwrite(uint8(I_temp(:,:,mm)), file_temp, 'jp2', 'CompressionRatio', new_compression_ratio);
            I_out_data(:, :, mm, kk) = double(imread(file_temp, 'jp2')) / 255 * dynamic_range;
        end
        delete(file_temp);
    end
end

fprintf('Quality index calculation:\n')
[MR_out,qi_label,MR_idx]= indexes_evaluation_RR(I_out_data,I_GT,'qindex_list',qindex_list);
            
output_folder=fullfile(output_folder,im_tag);
mkdir(output_folder);

filename=sprintf('%s',im_tag);

matrix2latex(MR_out.','filename',fullfile(output_folder,[filename,'.tex']),...
'row',methods_list,'col',qi_label,'align','c','significant',4,...
'bold',MR_idx.'==1,'underline',MR_idx.'==2); 

if flag_vis==1
    a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display));
    imwrite(a,fullfile(output_folder,[filename,'_GT.png']));
    for mm = 1:size(I_out_data,4)
        a=viewimage_outputonly(I_out_data(:, :, I_GT.Bands_to_display, mm));
        filename_current = sprintf([filename,'_%s.png'], methods_list{mm});
        imwrite(a,fullfile(output_folder, filename_current));
    end
end

I_out = I_GT;
I_out.data = I_out_data;
I_out.methods = methods_list;
I_out.label = 'SoftwareCompression';

MR.data = MR_out;
MR.qindex = qi_label;
MR.label = methods_list;
MR.bestindex = MR_idx;

save(fullfile(output_folder,[filename,'.mat']),'I_out','I_GT','MR');