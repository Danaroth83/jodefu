function PAN_LP = genPAN_LP_Equaliz(I_PAN,I_MS,equaliz, mode, ratio, varargin)
% Create a low-pass version of the PAN. The output can be a single or
% multi-band PAN at a lower resolution.
%
% Examples
%   PAN_LP = genPAN_LP(PAN, 'GLP', 4);
%
%   h = genMTF(4, 'QB');
%   PAN_LP = genPAN_LP(PAN, 'GLP_MTF', 4, h);
% keyboard
I_PAN = double(I_PAN);

if equaliz==0,
    imageHR = I_PAN;
else
    if size(I_PAN,3)==1
        I_PAN=repmat(I_PAN,[1 1 size(I_MS,3)]);
    end
    switch equaliz
        case 1
            for ii = 1:size(I_PAN,3),
                imageHR(:,:,ii) = I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii))+ mean2(I_MS(:,:,ii));
            end
        case 2
            for ii = 1:size(I_PAN,3),
                imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                    std2(I_MS(:,:,ii))/std2(I_PAN(:,:,ii)) + mean2(I_MS(:,:,ii));
            end
        case 3
            for ii = 1:size(I_PAN,3),
                imageHR(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii)))*...
                    std2(I_MS(:,:,ii))/std2(imresize(I_PAN(:,:,ii),1/ratio)) +...
                    mean2(I_MS(:,:,ii));
            end
    end
end

switch mode

    case 'GLP'
        fprintf('GLP [Aiazzi02]\n');
        
        PAN_LP = deg23tap(imageHR,ratio);
        PAN_LP = interp23tap(PAN_LP,ratio);
        
    case 'GLP_MTF_PAN'
        fprintf('GLP MTFPAN-matched [Aiazzi06]\n');
        
        if nargin < 4
            error('For mode %s there should be as an additional input: the filters windows', mode);
        end
        h = varargin{1};
        
        PAN_LP = repmat(imageHR, [1, 1, size(h,3)]);
        for i=1:size(h, 3)
            PAN_LP(:,:,i) = imfilter(imageHR,real(h(:,:,i)),'replicate');
        end
        
        %         PAN_LP = double(PAN_LP);
        PAN_LP = imresize(PAN_LP,1/ratio,'nearest');
        PAN_LP = interp23tap(PAN_LP,ratio);
    case 'GLP_MTF'
        fprintf('GLP MTF-matched [Aiazzi06]\n');
        
        if nargin < 6
            error('For mode %s there should be as an additional input: the filters windows', mode);
        end
        h = varargin{1};
        if size(imageHR,3)~=size(h,3)
            imageHR = repmat(imageHR, [1, 1, size(h,3)]);
        end
        for i=1:size(h, 3)
            PAN_LP(:,:,i) = imfilter(imageHR(:,:,i),real(h(:,:,i)),'replicate');
        end
        
        %         PAN_LP = double(PAN_LP);
        PAN_LP = imresize(PAN_LP,1/ratio,'nearest');
        if (2^round(log2(ratio)) ~= ratio)
            PAN_LP = interp23tapGeneral(PAN_LP,ratio);
        else
            PAN_LP = interp23tap(PAN_LP,ratio);
        end
    case 'undecimated_GLP_MTF'
%         keyboard
        fprintf('GLP MTF-matched [Aiazzi06]\n');
        
        if nargin < 4
            error('For mode %s there should be as an additional input: the filters windows', mode);
        end
        h = varargin{1};
        
        if size(imageHR,3)~=size(h,3)
            imageHR = repmat(imageHR, [1, 1, size(h,3)]);
        end
        
        for i=1:size(h, 3)
            PAN_LP(:,:,i) = imfilter(imageHR(:,:,i),real(h(:,:,i)),'replicate');
        end
        
        %         PAN_LP = double(PAN_LP);
        
    case 'ATWT'
        
        h=[1 4 6 4 1 ]/16;
        g=[0 0 1 0 0 ]-h;
        htilde=[ 1 4 6 4 1]/16;
        gtilde=[ 0 0 1 0 0 ]+htilde;
        h=sqrt(2)*h;
        g=sqrt(2)*g;
        htilde=sqrt(2)*htilde;
        gtilde=sqrt(2)*gtilde;
        WF={h,g,htilde,gtilde};
        
        Levels = ceil(log2(ratio));
        
        WT = ndwt2(imageHR,Levels,WF);
        
        for ii = 2 : numel(WT.dec), WT.dec{ii} = zeros(size(WT.dec{ii})); end
        
        PAN_LP = indwt2(WT,'c');
        
        %         PAN_LP = interp23tap(LPfilterPlusDec(I_PAN,ratio),ratio);
        
end

