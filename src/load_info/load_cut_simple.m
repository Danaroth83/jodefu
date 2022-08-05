function [vcut,hcut,edge_cut_lq,edge_cut_hq,Qblocks_size]=load_cut_simple(im_tag,cut_label)

if strncmpi(cut_label, 'cut', 3)
    vertical_length = str2double(cut_label(4:end));
else
    vertical_length = 512;
end
horizontal_length = vertical_length;

if strcmpi(im_tag, 'Hobart')
    edge_cut_lq = floor((512 - vertical_length) / 2);
else
    edge_cut_lq = floor((600 - vertical_length) / 2);
end
if edge_cut_lq < 0
   error('Requested cropped area is bigger than the image size');
end

vcut = edge_cut_lq + (1:vertical_length);
hcut = edge_cut_lq + (1:horizontal_length);

edge_cut_lq = 0;
edge_cut_hq = 0;

Qblocks_size = 32;
    
end