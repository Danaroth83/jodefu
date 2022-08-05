%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Visualization [3-2-1] of images with 3 bands by exploiting linear stretching and fixing the saturation. 
% 
% Interface:
%           ImageToView = viewimage(ImageToView,tol)
%
% Inputs:
%           ImageToView:    Image to view;
%           tol:            Saturation; Default values: [0.01 0.99] equal for all the three bands.
%
% Outputs:
%           ImageToView:    Image to view.
%           tlim:           Limits for imshow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ImageToView,tlim] = viewimage_outputonly(ImageToView,tol)

    if nargin<=1 || isempty(tol), tol=[0.01, 0.99]; end
    if numel(tol)==2, tol=[tol(1),tol(2)]; tol=[tol;tol;tol]; end

    ImageToView = double(ImageToView);
    if size(ImageToView,1)*size(ImageToView,2)>1200^2, error('Image is too big to show'); end
    L=size(ImageToView,3);
    if L<3, ImageToView=ImageToView(:,:,[1 1 1]); end

    if sum(tol(:,2)) <= 3
        [ImageToView,tlim] = linstretch(ImageToView,tol);
    else
        [N,M,~] = size(ImageToView);
        NM = N*M;
        for ii=1:3
            b = reshape(double(uint16(ImageToView(:,:,ii))),NM,1);
            b(b<tol(ii,1))=tol(ii,1);
            b(b>tol(ii,2))=tol(ii,2);
            b = (b-tol(ii,1))/(tol(ii,2)-tol(ii,1));
            ImageToView(:,:,ii) = reshape(b,N,M);
        end
    end
    
    % iptsetpref('ImshowBorder', 'tight')
    % figure,imshow(ImageToView,[])
    % iptsetpref('ImshowBorder', 'loose');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Linear Stretching. 
% 
% Interface:
%           ImageToView = linstretch(ImageToView,tol)
%
% Inputs:
%           ImageToView:    Image to stretch;
%           tol:            Ratio of tolerance;
%
% Outputs:
%           ImageToView:    Stretched image.
%           t:              Limit values as absolute pixel values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ImageToView,t] = linstretch(ImageToView,tol)

    [N,M,~] = size(ImageToView);
    NM = N*M;
    t=zeros(3,2);
    for ii=1:3
        % b = reshape(double(uint16(ImageToView(:,:,i))),NM,1);
        b = reshape(double((ImageToView(:,:,ii))),NM,1);
        % [hb,levelb] = hist(b,max(b)-min(b));
        [hb,levelb] = hist(b,10000);

        chb = cumsum(hb);
        % t(ii,1)=ceil(levelb(find(chb>NM*tol(ii,1), 1 )));
        % t(ii,2)=ceil(levelb(find(chb<NM*tol(ii,2), 1, 'last' )));
        t(ii,1)=(levelb(find(chb>NM*tol(ii,1), 1 )));
        t(ii,2)=(levelb(find(chb<NM*tol(ii,2), 1, 'last' )));
        b(b<t(ii,1))=t(ii,1);
        b(b>t(ii,2))=t(ii,2);
        b = (b-t(ii,1))/(t(ii,2)-t(ii,1));
        ImageToView(:,:,ii) = reshape(b,N,M);
    end

end