%% IMAGE VISUALIZATION
%
% Description:
%
% This function plots a series of images in a appropriately contrast
% stretched intensity range
%
% Interface:
%
% I_vis=Visualization_results(Image_cell,Field,Value);
% eg: MI=Visualization_results({MS,PAN},'filename','Pansharpening\output','printMAT',1);
%
% Inputs:
% Im_cell: cell containing a series of struct whose fields must be:
%        Bands_to_display: the indices of the three bands to display
%        methods_list: list of labels of the image
%        label: type of image (if above is not present)
%
% Other Fields:
%     output_file: Name of the image file to save (all results are saved under '..\..\output\')
%     printEPS: If 1, print results as EPS (default=0)
%     printTIFF: If 1, print results as TIFF (default=1)
%     printMAT: If 1, print results as MAT (default=1)
%     flag_draw: If 1, keep results on screen (default=1)
%     flag_cutbands: If 1, remove dim_cut pixels from edges for visualization
%     dim_cut: amount of pixels to remove from visualization

function I_vis_out=Visualization_image(Im_cell,varargin)

%% Definition of folders
current_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(current_folder,'..','..','data','output');

%% Definition of variables
output_file='output';
printEPS = 0;
printTIFF = 1;
printMAT = 0;
flag_draw = 1;
flag_cutbounds=0;
dim_cut=0;
flag_thvalues=0;

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if any(strcmpi(pname,{'output_file','output_name','filename'}))
        idx=strfind(pval,'.');
        if isempty(idx) || idx(end)<=1, output_file=pval; else, output_file=pval(1:idx(end)-1); end
    elseif any(strcmpi(pname,{'EPS','printEPS'}))   % If 1, print results as EPS
        printEPS=pval;
    elseif any(strcmpi(pname,{'TIFF','printTIFF'})) % If 1, print results as TIFF
        printTIFF=pval;
    elseif any(strcmpi(pname,{'flag_keep','draw','flag_draw'})) % Shift to mask
        flag_draw=pval;
    elseif any(strcmpi(pname,{'flag_cutbounds','flag_cutbounds_vis','cutbounds'})) % If 1, it removes edges in visualization
        flag_cutbounds=pval; 
    elseif any(strcmpi(pname,{'dim_cut','dim_cut_vis'})) % Shift to mask
        dim_cut=pval;
    elseif any(strcmpi(pname,{'flag_thvalues','flag_thvalues_vis','thvalues'})) % Shift to mask
        flag_thvalues=pval;
    end
end

output_fullname=fullfile(output_folder,output_file);

if ~iscell(Im_cell), Im_cell={Im_cell}; flag_imcell=1; else, flag_imcell=0; end

Xmap=10; Ymap=19;

N_cells=numel(Im_cell);
I_vis_out=cell(1,N_cells);

for jj=1:N_cells
    if isfield(Im_cell{jj},'image'), I_load=Im_cell{jj}.image; else, I_load=Im_cell{jj}; end
    if isfield(Im_cell{jj},'Bands_to_display'), Bands_to_display=Im_cell{jj}.Bands_to_display; else, Bands_to_display=1:min(3,size(I_load,3)); end
    if isfield(Im_cell{jj},'IntBits'), L=Im_cell{jj}.IntBits; else, L=8; end
    if isfield(Im_cell{jj},'DynamicRange'), DynamicRange=Im_cell{jj}.DynamicRange; else, DynamicRange=2^L-1; end
    if isfield(Im_cell{jj},'ratio'), ratio=Im_cell{jj}.ratio; else, ratio=1; end
    if isfield(Im_cell{jj},'titleImages')
        methods_list=Im_cell{jj}.titleImages;
    elseif isfield(Im_cell{jj},'methods_list')
        methods_list=Im_cell{jj}.methods_list;
    elseif isfield(Im_cell{jj},'label')
        methods_list={Im_cell{jj}.label};
    else
        for ii=1:size(Im_cell{jj}.image,4), methods_list{ii}=sprintf('UNK%d',ii); end
    end
    
    if flag_cutbounds
        I_load = I_load(round(dim_cut/ratio):end-round(dim_cut/ratio),round(dim_cut/ratio):end-round(dim_cut/ratio),:,:);
    end

    if flag_thvalues
        I_load(I_load > DynamicRange) = DynamicRange;
        I_load(I_load < 0) = 0;
    end
    
    [L1,L2,~,Nt]=size(I_load);
    I_out_temp=zeros([L1,L2,3,Nt]);
    
    for kk=1:Nt
        
        method_filename=strrep(methods_list{kk},'-','_');        
        I_vis=I_load(:,:,Bands_to_display,kk);
        
        
        I_out_temp(:,:,:,kk)=viewimage_outputonly(I_vis);
        
        if printTIFF==1 || printEPS==1
            figure;
            imshow(I_out_temp(:,:,:,kk) ,'Border','tight','InitialMagnification',100);
        end
        
        if printTIFF==1
            filename_TIFF=[output_fullname,'_',method_filename,'.tif'];
            option_plot=get(gcf,'PaperPositionMode');
            set(gcf,'PaperPositionMode','auto');
            print(filename_TIFF,'-dtiffn','-r0');
            set(gcf,'PaperPositionMode',option_plot);
        end
        if printEPS==1
            filename_EPS=[output_fullname,'_',method_filename,'.eps'];
            print(filename_EPS,'-depsc2','-r300');
        end
        if flag_draw==1
            method_label=strrep(methods_list{kk},'_','-'); 
            text(Xmap,Ymap,method_label,'EdgeColor','k','BackgroundColor','w');
        else
            close;
        end
    end
    I_vis_out{jj}=I_out_temp;
end

if flag_imcell==0, I_vis_out=I_vis_out{1}; end

if printMAT==1
    save([output_fullname,'.mat'],'I_vis_out','-append');
end

drawnow;