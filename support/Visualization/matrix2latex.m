function textstring=matrix2latex(matrix, varargin)

% function: matrix2latex(...)
% Author:   M. Koehler
% Contact:  koehler@in.tum.de
% Version:  1.4 (Modified by Daniele Picone)
% Date:     April 12, 2020

% This software is published under the GNU GPL, by the free software
% foundation. For further reading see: http://www.gnu.org/licenses/licenses.html#GPL

% Usage:
% matrix2latex(matrix varargs)
% where
%   - matrix is a 2 dimensional numerical or cell array
%   
%   - varargs is one ore more of the following (denominator, value) combinations
%   + filename is a valid filename, in which the resulting latex code will
%      be stored (within the directory ..\..\output\')
%      + 'rowLabels', array -> Can be used to label the rows of the
%      resulting latex table
%      + 'columnLabels', array -> Can be used to label the columns of the
%      resulting latex table
%      + 'alignment', 'value' -> Can be used to specify the alginment of
%      the table within the latex document. Valid arguments are: 'l', 'c',
%      and 'r' for left, center, and right, respectively
%      + 'format', 'value' -> Can be used to format the input data. 'value'
%      has to be a valid format string, similar to the ones used in
%      fprintf('format', value);
%      + 'size', 'value' -> One of latex' recognized font-sizes, e.g. tiny,
%      HUGE, Large, large, LARGE, etc.
%      + 'highlight', 'value': highliting of the best results ('n': no
%      highlighting, 'g': higlighting highest values, 's': smallest)
%      + 'gRow', grouping
%
% Example input:
%   matrix = [1.5 1.764; 3.523 0.2];
%   rowLabels = {'row 1', 'row 2'};
%   columnLabels = {'col 1', 'col 2'};
%   matrix2latex(matrix, 'output', 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%
% The resulting latex file can be included into any latex document by:
% /input{out.tex}
%
% Enjoy life!!!
    
width = size(matrix, 2);
height = size(matrix, 1);

rowLabels = [];
colLabels = [];
alignment = 'c';
format = [];
textsize = [];
highlight='n';
groups_r=height;
groups_c=width;
groups_hr=height;
groups_hc=ones(1,width);
significant=4;
filename=fullfile(fileparts(mfilename('fullpath')),'..','..','output','output.tex');
printTEX=0;
highlight_matrix_bold=zeros(size(matrix));
highlight_matrix_underline=zeros(size(matrix));
flag_highlight=0;

if (rem(nargin,2) == 0 || nargin < 1)
    error(['matrix2latex: ', 'Incorrect number of arguments to %s.'], mfilename);
end

for j=1:2:(nargin-2)
    pname = varargin{j};
    pval = varargin{j+1};
    if any(strcmpi(pname,{'rowlabels','row'}))  % rowlabels
        rowLabels = pval;
        if isnumeric(rowLabels)
            rowLabels = cellstr(num2str(rowLabels(:)));
        end
    elseif any(strcmpi(pname,{'columnlabels','column','col'}))  % column labels
        colLabels = pval;
        if isnumeric(colLabels)
            colLabels = cellstr(num2str(colLabels(:)));
        end
    elseif any(strcmpi(pname,{'alignment','align'}))  % alignment
        alignment = lower(pval);
        if strcmp(alignment,'right')
            alignment = 'r';
        end
        if strcmp(alignment,'left')
            alignment = 'l';
        end
        if strcmp(alignment,'center')
            alignment = 'c';
        end
        if alignment ~= 'l' && alignment ~= 'c' && alignment ~= 'r'
            alignment = 'c';
            warning(['matrix2latex: ', 'Unknown alignment. (Set it to \''center\''.)']);
        end
    elseif strcmpi(pname,{'format'})  % format
        format = lower(pval);
    elseif any(strcmpi(pname,{'size','textsize'}))  % text size
        textsize = pval;
    elseif strcmpi(pname,{'highlight'})  % highlight
        highlight=pval;
        if highlight(1)~='c' && highlight(1)~='r' && highlight(1)~='n'
            error('First character must be ''c'' for highlighting columns, ''r'' for rows, ''n'' for no alignment');
        else
            strrep(highlight(2:end),'l','s');
            strrep(highlight(2:end),'b','g');
            if ~all(highlight(2:end)=='s' | highlight(2:end)=='g')
                error('Use ''g'' for highlighting greatest value, ''s'' for smallest');
            end
        end
    elseif any(strcmpi(pname,{'gRow','groups_r'})) % Grouping of rows
        if sum(pval)~=height, error('Row grouping is wrong'); end
        groups_r=ismember(1:height,cumsum(pval));
    elseif any(strcmpi(pname,{'gCol','groups_c'})) % Grouping of columns
        if sum(pval)~=width, error('Column grouping is wrong'); end
        groups_c=ismember(1:width,cumsum(pval));
    elseif any(strcmpi(pname,{'hRow','groups_hr'})) % Grouping of highlights for rows
        if sum(pval)~=height, error('Row highlighting grouping is wrong'); end
        groups_hr=pval;
    elseif any(strcmpi(pname,{'hCol','groups_hc'})) % Grouping of highlights for columns
        if sum(pval)~=width, error('Column highlighting grouping is wrong'); end
        groups_hc=pval;
    elseif any(strcmpi(pname,{'flag_highlight'})) % Standard ways of higlighting
        flag_highlight=pval;
    elseif any(strcmpi(pname,{'significant','digits'})) % Significant digits
        significant=pval;
    elseif any(strcmpi(pname,{'highlight_bold','bold'}))
        highlight_matrix_bold=pval;
    elseif any(strcmpi(pname,{'highlight_underline','underline'}))
        highlight_matrix_underline=pval;    
    elseif any(strcmpi(pname,{'filename','file','name','output','output_name'})) % filename
        idx=strfind(pval,'.');
        if isempty(idx) || idx(end)<=1, filename=pval; else, filename=pval(1:idx(end)-1); end
        printTEX=1;
    elseif any(strcmpi(pname,{'printTEX','TEX'}))
        printTEX=pval;
    else
        error('Unknown Parameter');
    end
end

if flag_highlight==0
    highlight='n';
elseif flag_highlight==1
    groups_hc=ones(1,width);
    if isempty(colLabels), error('Field colLabels is required'); end
    highlight=matrix2latex_highlightoption(colLabels);   
end

groups_c=cumsum(groups_c); groups_c_temp=zeros(1,width); groups_c_temp(groups_c)=1; groups_c_temp(end)=0;
groups_r=cumsum(groups_r); groups_r_temp=zeros(1,height); groups_r_temp(groups_r)=1; groups_r_temp(end)=0;
groups_hc=cumsum(groups_hc);
groups_hr=cumsum(groups_hr);

if (highlight(1)=='c' && groups_hc(end)~=width) || (highlight(1)=='r' && groups_hr(end)~=height)
    error('Highlight grouping is wrong');
end
if (highlight(1)=='c' && length(groups_hc)~=length(highlight)-1) || (highlight(1)=='r' && length(groups_hr)~=length(highlight)-1)
    error('Highlighting dimensions are not correct');
end

% matrix=round(matrix,significant,'significant');
% matrix_orig=matrix;
% if ~isempty(significant), round(matrix,significant,'significant'); end

if isnumeric(matrix)
    % matrix = num2cell(matrix);
    matrix_num=matrix;
    matrix=cell(size(matrix));
    for h=1:height
        for w=1:width
            if (~isempty(significant))
                format=['%.',sprintf('%d',max(-floor(log10(matrix_num(h,w)))+significant-1,0)),'f'];
                matrix{h, w} = num2str(matrix_num(h,w),format);
            elseif (~isempty(format))
                matrix{h, w} = num2str(matrix_num(h,w),format);
            else
                matrix{h, w} = num2str(matrix_num(h,w));
            end
        end
    end
end


if ~strcmpi(highlight,'n')
    start_groups_hr=[1,groups_hr+1]; end_groups_hr=start_groups_hr(2:end)-1;
    start_groups_hc=[1,groups_hc+1]; end_groups_hc=start_groups_hc(2:end)-1;
    for i1=1:length(groups_hr)
        idx1=start_groups_hr(i1):end_groups_hr(i1);
        for i2=1:length(groups_hc)
            idx2=start_groups_hc(i2):end_groups_hc(i2);
            matrix_block=matrix_orig(idx1,idx2);
            if highlight(1)=='r'
                ordering=highlight(i1+1);
            else
                ordering=highlight(i2+1);
            end
            if ordering=='g'
                matrix_ordered=sort(matrix_block(:),'descend');
            else
                matrix_ordered=sort(matrix_block(:),'ascend');
            end
            if numel(matrix_block)>1
                highlight_matrix_bold(idx1,idx2)=(matrix_block==matrix_ordered(1));
            end
            if numel(matrix_block)>2
                highlight_matrix_underline(idx1,idx2)=(matrix_block==matrix_ordered(2));
            end
        end
    end
end


textstring=[];

if(~isempty(textsize))
    textstring=cat(2,textstring, ['\\begin{',sprintf('%s', textsize),'}']);
end

textstring=cat(2,textstring,'\\begin{tabular}{');

if(~isempty(rowLabels))
    textstring=cat(2,textstring,'l|');
    textstring=cat(2,textstring,sprintf('%s',repmat('l|',[1,length(strfind(colLabels{1},'&'))])));
end
for w=1:width
        textstring=cat(2,textstring,sprintf('%c', alignment));
    if groups_c_temp(w)==1
        textstring=cat(2,textstring,'|');
    end
end
textstring=cat(2,textstring,'}\r\n');

if(~isempty(colLabels))
    if(~isempty(rowLabels))
        textstring=cat(2,textstring,'&');
    end
    for w=1:width
        title=colLabels{w};
        idx=strfind(title,'&');
        if ~isempty(idx), idx=idx(end); else, idx=0; end
        if isempty(regexp(title(idx+1:end),'[_^=]','once'))
            textstring=cat(2,textstring,[title(1:idx),'\\textbf{',sprintf('%s',title(idx+1:end)),'}']);
        else
            textstring=cat(2,textstring,[title(1:idx),'$\\mathbf{',sprintf('%s',title(idx+1:end)),'}$']);
        end
        if w==width
            textstring=cat(2,textstring,'\\\\\\hline\r\n');
        else
            textstring=cat(2,textstring,'&');
        end
    end
end


for h=1:height
    if(~isempty(rowLabels))
        title=rowLabels{h};
        idx=strfind(title,'&');
        if ~isempty(idx), idx=idx(end); else, idx=0; end
        if isempty(regexp(title(idx+1:end),'[_^=]','once'))
            textstring=cat(2,textstring,[title(1:idx),'\\textbf{',sprintf('%s',title(idx+1:end)),'}&']);
        else
            textstring=cat(2,textstring,[title(1:idx),'$\\mathbf{',sprintf('%s',title(idx+1:end)),'}$&']);
        end
    end
    for w=1:width
        if highlight_matrix_bold(h,w)==1
            textstring=cat(2,textstring,['\\textbf{',sprintf('%s',matrix{h, w}),'}']);
        elseif highlight_matrix_underline(h,w)==1
            textstring=cat(2,textstring,['\\underline{',sprintf('%s',matrix{h, w}),'}']);
        else
            textstring=cat(2,textstring,sprintf('%s', matrix{h, w}));
        end
        if w==width
            if groups_r_temp(h)==1
                textstring=cat(2,textstring,'\\\\\\hline\r\n');
            else
                textstring=cat(2,textstring,'\\\\\r\n');
            end
        else
            textstring=cat(2,textstring,'&');
        end
    end
end

textstring=cat(2,textstring,'\\end{tabular}\r\n');

if(~isempty(textsize))
    textstring=cat(2,textstring,['\\end{',sprintf('%s', textsize),'}']);
end

if printTEX==1
    fid = fopen(fullfile([filename,'.tex']), 'a');
    fprintf(fid,textstring);
    fclose(fid);
end