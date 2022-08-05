function output = iswt2_Bspline(a,ratio,order)
%ISWT2 Inverse discrete stationary wavelet transform 2-D.
%   ISWT2 performs a multilevel 2-D stationary wavelet 
%   reconstruction using either a specific orthogonal wavelet  
%   ('wname', see WFILTERS for more information) or specific
%   reconstruction filters (Lo_R and Hi_R).
%
%   X = ISWT2(SWC,'wname') or X = ISWT2(A,H,V,D,'wname') 
%   or X = ISWT2(A(:,:,end),H,V,D,'wname') reconstructs the 
%   matrix X, based on the multilevel stationary wavelet   
%   decomposition structure SWC or [A,H,V,D] (see SWT2).
%
%   For X = ISWT2(SWC,Lo_R,Hi_R) or X = ISWT2(A,H,V,D,Lo_R,Hi_R) 
%   or X = ISWT2(A(:,:,end),H,V,D,Lo_R,Hi_R): 
%   Lo_R is the reconstruction low-pass filter.
%   Hi_R is the reconstruction high-pass filter.
%
%   NOTE: If  SWC or (CA, CH, CV, CD) are obtained from an 
%   indexed image analysis (respectively a truecolor image
%   analysis) then X is an m-by-n matrix (respectively an
%   m-by-n-by-3 array).
%   For more information on image formats, see the reference
%   pages of IMAGE and IMFINFO functions.
%   
%   See also IDWT2, SWT2, WAVEREC2.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 08-Dec-97.
%   Last Revision: 20-Dec-2010.
%   Copyright 1995-2010 The MathWorks, Inc.

% Check arguments.
% nbIn = nargin;
% switch nbIn
%     case {0,1}
%         error(message('Wavelet:FunctionInput:NotEnough_ArgNum'));
%     case 4
%         error(message('Wavelet:FunctionInput:Invalid_ArgNum'));
%     case {2,3,5,6}
%     otherwise
%         error(message('Wavelet:FunctionInput:TooMany_ArgNum'));
% end
% switch nbIn
%     case 2 , argstr = 1; argnum = 2;
%     case 3 , argstr = 0; argnum = 2;
%     case 5 , argstr = 1; argnum = 5;
%     case 6 , argstr = 0; argnum = 5;
% end
% 
% % Compute decomposition filters.
% if argstr
%     [lo_R,hi_R] = wfilters(varargin{argnum},'r');
% else
%     lo_R = varargin{argnum}; hi_R = varargin{argnum+1};
% end

if nargin <=2
    order=3;
end

% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

% % Get inputs.
% sX = size(a);
% level=sX(3);
% lenSX = length(sX);
% a3d_Flag = lenSX>3 || ((lenSX==3) && isequal(class(a),'single'));
% if a3d_Flag , nbColon = 3; else nbColon = 2; end  % max([lenSX-1,2])
% idxColon = repmat({':'},1,nbColon);
% S.type = '()';
% Ambiguity = false;
% if argnum==2
%     nbRows = size(a,lenSX);
%     if rem(nbRows,3)==1
%         level = (nbRows-1)/3;
%     else
%         errargt(mfilename,'Invalid size for the first argument!','msg');
%         error(message('Wavelet:FunctionArgVal:Invalid_Input'));
%     end
%     S.subs = [idxColon,1:level];
%     h = subsref(varargin{1},S);
%     S.subs = [idxColon,level+1:2*level];
%     v = subsref(varargin{1},S);
%     S.subs = [idxColon,2*level+1:3*level];
%     d = subsref(varargin{1},S);
%     S.subs = [idxColon,3*level+1:nbRows];
%     a = subsref(varargin{1},S);
% else
%     a = varargin{1};
%     h = varargin{2};
%     v = varargin{3};
%     d = varargin{4};
%     %----------------------------------------------------------------------
%     % Ambiguity:
%     % (level=3 and indexed BW image) or (level=1 and truecolor image)
%     % To suppress this Ambiguity, the function SWT2 in case of a true 
%     % color image and a level 1 analysis, produce single for approximation
%     % coefficients !!
%     %----------------------------------------------------------------------
%     % if ~a3d_Flag
%     %     if size(a,3)==3 && size(h,3)==3 && size(v,3)==3 && size(d,3)==3
%     %         a3d_Flag = true;
%     %         Ambiguity = true;
%     %     end
%     % end
%     %----------------------------------------------------------------------    
% end
[rx,cx,dim3,dim4] = size(a);

% Extract last approximation coefficients
% if ~Ambiguity
%     if a3d_Flag
%         idxApp = size(a,4);
%     else
%         idxApp = size(a,3);
%     end
%     S.subs = [idxColon,idxApp];
%     a  = subsref(a,S);
% else  % do nothing    
%     error(message('Wavelet:FunctionToVerify:LastApp'));
% end

% if ~a3d_Flag
%     level = dim3;
%     a = reconsLOC(a,h,v,d);    
% else
%     level = dim4;
%     tmp = cell(1,3);
%     for j=1:3         
%         tmp{j} = reconsLOC(a(:,:,j),h(:,:,j,:),v(:,:,j,:),d(:,:,j,:));
%     end
%     a = cat(3,tmp{:});
%     a(a<0) = 0;
%     a = uint8(a);   
% end

a=reconsLOC(a,ratio,order);
output = a;

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');

%--------------------------------------------------------------------
    function ca = reconsLOC(ca,ratio,order)
        fac=factor(ratio);
        fac_prod=[1,cumprod(fac)];
        level=length(fac);
        for k = level:-1:1
            fac_k=fac(k);
            step = fac_prod(k);
            last = step;
            
            st_ones=ones(1,fac_k);
            lo_R=st_ones;
            for ii=1:order
                lo_R=conv(lo_R,st_ones);
            end
            lo_R=sqrt(fac_k)/fac_k^(order+1)*lo_R;

            % %upsample filter
            % lo_R(2:fac_prod(k),length(lo_R))=0;
            % lo_R=lo_R(:)';
            
            for first = 1:last
                iRow = first:step:rx;
                lR   = length(iRow);
                iCol = first:step:cx;
                lC   = length(iCol);
                
                x=ones(lR,lC);
                for ii=1:fac_k
                    sR = iRow(ii:fac_k:lR);
                    sC = iCol(ii:fac_k:lC);
                    x = x + idwt2LOC(ca(sR,sC),lo_R,[lR,lC],[-(ii-1),-(ii-1)],fac_k);
                end
                ca(iRow,iCol)=x/fac_k;
%                 sR   = iRow(1:2:lR);
%                 sC   = iCol(1:2:lC);
%                 x1   = idwt2LOC(...
%                             ca(sR,sC),ch(sR,sC,k),cv(sR,sC,k),cd(sR,sC,k), ...
%                             lo_R,hi_R,[lR lC],[0,0]);
% 
%                 sR   = iRow(2:2:lR);
%                 sC   = iCol(2:2:lC);
%                 x2   = idwt2LOC(...
%                             ca(sR,sC),ch(sR,sC,k),cv(sR,sC,k),cd(sR,sC,k), ...
%                             lo_R,hi_R,[lR lC],[-1,-1]);
%                 ca(iRow,iCol) = 0.5*(x1+x2);
            end
        end
    end
%--------------------------------------------------------------------

end


%===============================================================%
% INTERNAL FUNCTIONS
%===============================================================%
function y = idwt2LOC(a,lo_R,sy,shift,curr_ratio)

y = upconvLOC(a,lo_R,lo_R,sy,curr_ratio); % Approximation.

y = y([end+shift(1)+1:end,1:end+shift(1)],[end+shift(2)+1:end,1:end+shift(2)]);
end
%---------------------------------------------------------------%
function y = upconvLOC(x,f1,f2,s,curr_ratio)

lf = length(f1);
% y  = dyadup(x,'mat',0,1);
y=zeros(curr_ratio*size(x));
y(1:curr_ratio:end,1:curr_ratio:end)=x;
y  = wextend('2D','per',y,[floor(lf/2),floor(lf/2)]);
y  = wconv2('col',y,f1);
y  = wconv2('row',y,f2);
y  = wkeep2(y,s,[lf lf]);
end
%===============================================================%
