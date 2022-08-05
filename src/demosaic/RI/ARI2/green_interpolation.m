function green = green_interpolation(mosaic, mask, pattern, eps, DynamicRange)

% inverse mask
imask = (mask == 0);

% raw CFA data
rawq = sum(mosaic, 3);

% maskGr: mask of G pixels horizontally neighboring R
% maskGb: mask of G pixels horizontally neighboring B
maskGr = zeros(size(rawq));
maskGb = zeros(size(rawq));
if strcmp(pattern,'grbg')
    maskGr(1:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
    maskGb(2:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'rggb')
    maskGr(1:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
    maskGb(2:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'gbrg')
    maskGb(1:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
    maskGr(2:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
elseif strcmp(pattern,'bggr')
    maskGb(1:2:size(rawq,1), 2:2:size(rawq,2)) = 1;
    maskGr(2:2:size(rawq,1), 1:2:size(rawq,2)) = 1;
end
% Mrh: mask of horizontal R and G line
% Mbh: mask of horizontal B and G line
% Mrv: mask of vertical R and G line
% Mbv: mask of vertical B and G line
Mrh = mask(:,:,1) + maskGr;
Mbh = mask(:,:,3) + maskGb;
Mrv = mask(:,:,1) + maskGb;
Mbv = mask(:,:,3) + maskGr;


%% Step (i): iterative directional interpolation
% initial linear interpolation (Eq.(5))
Kh = [1/2,0,1/2];
Kv = Kh';
rawh = imfilter(rawq, Kh, 'replicate');
rawv = imfilter(rawq, Kv, 'replicate');
% horizontal direction
Guidegrh = mosaic(:,:,2).*maskGr + rawh.*mask(:,:,1);
Guidegbh = mosaic(:,:,2).*maskGb + rawh.*mask(:,:,3);
Guiderh  = mosaic(:,:,1) + rawh.*maskGr;
Guidebh  = mosaic(:,:,3) + rawh.*maskGb;
% vertical direction
Guidegrv = mosaic(:,:,2).*maskGb + rawv.*mask(:,:,1);
Guidegbv = mosaic(:,:,2).*maskGr + rawv.*mask(:,:,3);
Guiderv  = mosaic(:,:,1) + rawv.*maskGb;
Guidebv  = mosaic(:,:,3) + rawv.*maskGr;

% initial guided filter window size for RI
h = 2;
v = 1;
% initial guided filter window size for MLRI
h2 = 4;
v2 = 0;
% maximum iteration number
itnum = 11;

% initialization of horizontal and vertical iteration criteria
  RI_w2h = ones(size(maskGr)) * 1e32;
  RI_w2v = ones(size(maskGr)) * 1e32;
MLRI_w2h = ones(size(maskGr)) * 1e32;
MLRI_w2v = ones(size(maskGr)) * 1e32;

% initial guide image for RI
  RI_Guidegrh = Guidegrh;
  RI_Guidegbh = Guidegbh;
  RI_Guiderh  = Guiderh;
  RI_Guidebh  = Guidebh;
  RI_Guidegrv = Guidegrv;
  RI_Guidegbv = Guidegbv;
  RI_Guiderv  = Guiderv;
  RI_Guidebv  = Guidebv; 
% initial guide image for MLRI
MLRI_Guidegrh = Guidegrh;
MLRI_Guidegbh = Guidegbh;
MLRI_Guiderh  = Guiderh;
MLRI_Guidebh  = Guidebh;
MLRI_Guidegrv = Guidegrv;
MLRI_Guidegbv = Guidegbv;
MLRI_Guiderv  = Guiderv;
MLRI_Guidebv  = Guidebv;

% initialization of interpolated G values
  RI_Gh = Guidegrh + Guidegbh;
  RI_Gv = Guidegrv + Guidegbv;
MLRI_Gh = Guidegrh + Guidegbh;
MLRI_Gv = Guidegrv + Guidegbv;

%% Iterative horizontal and vertical interpolation
for ittime = 1:itnum	
	% generate horizontal and vertical tentative estimate by RI (Eq.(6))
	RI_tentativeGrh = guidedfilter(RI_Guiderh, RI_Guidegrh, Mrh, h, v, eps);
	RI_tentativeGbh = guidedfilter(RI_Guidebh, RI_Guidegbh, Mbh, h, v, eps);
	RI_tentativeRh  = guidedfilter(RI_Guidegrh, RI_Guiderh, Mrh, h, v, eps);
	RI_tentativeBh  = guidedfilter(RI_Guidegbh, RI_Guidebh, Mbh, h, v, eps);
	RI_tentativeGrv = guidedfilter(RI_Guiderv, RI_Guidegrv, Mrv, v, h, eps);
	RI_tentativeGbv = guidedfilter(RI_Guidebv, RI_Guidegbv, Mbv, v, h, eps);
	RI_tentativeRv  = guidedfilter(RI_Guidegrv, RI_Guiderv, Mrv, v, h, eps);
	RI_tentativeBv  = guidedfilter(RI_Guidegbv, RI_Guidebv, Mbv, v, h, eps);
	% generate horizontal tentative estimate by MLRI
	Fh = [-1,0,2,0,-1];
	difR  = imfilter(MLRI_Guiderh ,Fh,'replicate');
	difGr = imfilter(MLRI_Guidegrh,Fh,'replicate');
	difB  = imfilter(MLRI_Guidebh, Fh,'replicate');
	difGb = imfilter(MLRI_Guidegbh,Fh,'replicate');
	MLRI_tentativeRh  = guidedfilter_MLRI(MLRI_Guidegrh, MLRI_Guiderh, Mrh, difGr, difR, mask(:,:,1), h2, v2, eps);
	MLRI_tentativeBh  = guidedfilter_MLRI(MLRI_Guidegbh, MLRI_Guidebh, Mbh, difGb, difB, mask(:,:,3), h2, v2, eps);
	MLRI_tentativeGrh = guidedfilter_MLRI(MLRI_Guiderh, MLRI_Guidegrh, Mrh, difR, difGr, maskGr     , h2, v2, eps);
	MLRI_tentativeGbh = guidedfilter_MLRI(MLRI_Guidebh, MLRI_Guidegbh, Mbh, difB, difGb, maskGb     , h2, v2, eps);
	% generate vertical tentative estimate by MLRI
	Fv = [-1;0;2;0;-1];
	difR  = imfilter(MLRI_Guiderv ,Fv,'replicate');
	difGr = imfilter(MLRI_Guidegrv,Fv,'replicate');
	difB  = imfilter(MLRI_Guidebv ,Fv,'replicate');
	difGb = imfilter(MLRI_Guidegbv,Fv,'replicate');
	MLRI_tentativeRv  = guidedfilter_MLRI(MLRI_Guidegrv, MLRI_Guiderv, Mrv, difGr, difR, mask(:,:,1), v2, h2, eps);
	MLRI_tentativeBv  = guidedfilter_MLRI(MLRI_Guidegbv, MLRI_Guidebv, Mbv, difGb, difB, mask(:,:,3), v2, h2, eps);
	MLRI_tentativeGrv = guidedfilter_MLRI(MLRI_Guiderv, MLRI_Guidegrv, Mrv, difR, difGr, maskGb     , v2, h2, eps);
	MLRI_tentativeGbv = guidedfilter_MLRI(MLRI_Guidebv, MLRI_Guidegbv, Mbv, difB, difGb, maskGr     , v2, h2, eps);
	
	% calculate residuals of RI and MLRI (Eq.(9))
	  RI_residualGrh = ( mosaic(:,:,2) -   RI_tentativeGrh ) .* maskGr;
	  RI_residualGbh = ( mosaic(:,:,2) -   RI_tentativeGbh ) .* maskGb;
	  RI_residualRh  = ( mosaic(:,:,1) -   RI_tentativeRh  ) .* mask(:,:,1);
	  RI_residualBh  = ( mosaic(:,:,3) -   RI_tentativeBh  ) .* mask(:,:,3);
	  RI_residualGrv = ( mosaic(:,:,2) -   RI_tentativeGrv ) .* maskGb;
	  RI_residualGbv = ( mosaic(:,:,2) -   RI_tentativeGbv ) .* maskGr;
	  RI_residualRv  = ( mosaic(:,:,1) -   RI_tentativeRv  ) .* mask(:,:,1);
	  RI_residualBv  = ( mosaic(:,:,3) -   RI_tentativeBv  ) .* mask(:,:,3);
	MLRI_residualGrh = ( mosaic(:,:,2) - MLRI_tentativeGrh ) .* maskGr;
	MLRI_residualGbh = ( mosaic(:,:,2) - MLRI_tentativeGbh ) .* maskGb;
	MLRI_residualRh  = ( mosaic(:,:,1) - MLRI_tentativeRh  ) .* mask(:,:,1);
	MLRI_residualBh  = ( mosaic(:,:,3) - MLRI_tentativeBh  ) .* mask(:,:,3);
	MLRI_residualGrv = ( mosaic(:,:,2) - MLRI_tentativeGrv ) .* maskGb;
	MLRI_residualGbv = ( mosaic(:,:,2) - MLRI_tentativeGbv ) .* maskGr;
	MLRI_residualRv  = ( mosaic(:,:,1) - MLRI_tentativeRv  ) .* mask(:,:,1);
	MLRI_residualBv  = ( mosaic(:,:,3) - MLRI_tentativeBv  ) .* mask(:,:,3);

	% horizontal and vertical linear interpolation of residuals (Eq.(10))
	Kh = [1/2,1,1/2];
	RI_residualGrh   = imfilter(  RI_residualGrh, Kh, 'replicate');
	RI_residualGbh   = imfilter(  RI_residualGbh, Kh, 'replicate');
	RI_residualRh    = imfilter(  RI_residualRh , Kh, 'replicate');
	RI_residualBh    = imfilter(  RI_residualBh , Kh, 'replicate');
	MLRI_residualGrh = imfilter(MLRI_residualGrh, Kh, 'replicate');
	MLRI_residualGbh = imfilter(MLRI_residualGbh, Kh, 'replicate');
	MLRI_residualRh  = imfilter(MLRI_residualRh , Kh, 'replicate');
	MLRI_residualBh  = imfilter(MLRI_residualBh , Kh, 'replicate');
	Kv = Kh';
	RI_residualGrv   = imfilter(  RI_residualGrv, Kv, 'replicate');
	RI_residualGbv   = imfilter(  RI_residualGbv, Kv, 'replicate');
	RI_residualRv    = imfilter(  RI_residualRv , Kv, 'replicate');
	RI_residualBv    = imfilter(  RI_residualBv , Kv, 'replicate');
	MLRI_residualGrv = imfilter(MLRI_residualGrv, Kv, 'replicate');
	MLRI_residualGbv = imfilter(MLRI_residualGbv, Kv, 'replicate');
	MLRI_residualRv  = imfilter(MLRI_residualRv , Kv, 'replicate');
	MLRI_residualBv  = imfilter(MLRI_residualBv , Kv, 'replicate');

	% add tentative estimate (Eq.(11))
	  RI_Grh = (   RI_tentativeGrh +   RI_residualGrh ) .* mask(:,:,1);
	  RI_Gbh = (   RI_tentativeGbh +   RI_residualGbh ) .* mask(:,:,3);
	  RI_Rh  = (   RI_tentativeRh  +   RI_residualRh  ) .* maskGr;
	  RI_Bh  = (   RI_tentativeBh  +   RI_residualBh  ) .* maskGb;
	  RI_Grv = (   RI_tentativeGrv +   RI_residualGrv ) .* mask(:,:,1);
	  RI_Gbv = (   RI_tentativeGbv +   RI_residualGbv ) .* mask(:,:,3);
	  RI_Rv  = (   RI_tentativeRv  +   RI_residualRv  ) .* maskGb;
	  RI_Bv  = (   RI_tentativeBv  +   RI_residualBv  ) .* maskGr;
	MLRI_Grh = ( MLRI_tentativeGrh + MLRI_residualGrh ) .* mask(:,:,1);
	MLRI_Gbh = ( MLRI_tentativeGbh + MLRI_residualGbh ) .* mask(:,:,3);
	MLRI_Rh  = ( MLRI_tentativeRh  + MLRI_residualRh  ) .* maskGr;
	MLRI_Bh  = ( MLRI_tentativeBh  + MLRI_residualBh  ) .* maskGb;
	MLRI_Grv = ( MLRI_tentativeGrv + MLRI_residualGrv ) .* mask(:,:,1);
	MLRI_Gbv = ( MLRI_tentativeGbv + MLRI_residualGbv ) .* mask(:,:,3);
	MLRI_Rv  = ( MLRI_tentativeRv  + MLRI_residualRv  ) .* maskGr;
	MLRI_Bv  = ( MLRI_tentativeBv  + MLRI_residualBv  ) .* maskGr;

    %% Step(ii): adaptive selection of iteration at each pixel
    % calculate iteration criteria (Eq.(12))
	  RI_criGrh = (   RI_Guidegrh -   RI_tentativeGrh ) .* Mrh;
	  RI_criGbh = (   RI_Guidegbh -   RI_tentativeGbh ) .* Mbh;
	  RI_criRh  = (   RI_Guiderh  -   RI_tentativeRh  ) .* Mrh;
	  RI_criBh  = (   RI_Guidebh  -   RI_tentativeBh  ) .* Mbh;
	  RI_criGrv = (   RI_Guidegrv -   RI_tentativeGrv ) .* Mrv;
	  RI_criGbv = (   RI_Guidegbv -   RI_tentativeGbv ) .* Mbv;
	  RI_criRv  = (   RI_Guiderv  -   RI_tentativeRv  ) .* Mrv;
	  RI_criBv  = (   RI_Guidebv  -   RI_tentativeBv  ) .* Mbv;
	MLRI_criGrh = ( MLRI_Guidegrh - MLRI_tentativeGrh ) .* Mrh;
	MLRI_criGbh = ( MLRI_Guidegbh - MLRI_tentativeGbh ) .* Mbh;
	MLRI_criRh  = ( MLRI_Guiderh  - MLRI_tentativeRh  ) .* Mrh;
	MLRI_criBh  = ( MLRI_Guidebh  - MLRI_tentativeBh  ) .* Mbh;
	MLRI_criGrv = ( MLRI_Guidegrv - MLRI_tentativeGrv ) .* Mrv;
	MLRI_criGbv = ( MLRI_Guidegbv - MLRI_tentativeGbv ) .* Mbv;
	MLRI_criRv  = ( MLRI_Guiderv  - MLRI_tentativeRv  ) .* Mrv;
	MLRI_criBv  = ( MLRI_Guidebv  - MLRI_tentativeBv  ) .* Mbv;

    % calculate gradient of iteration criteria
	Fh = [-1,0,1];
	  RI_difcriGrh = abs(imfilter(  RI_criGrh, Fh, 'replicate'));
	  RI_difcriGbh = abs(imfilter(  RI_criGbh, Fh, 'replicate'));
	  RI_difcriRh  = abs(imfilter(  RI_criRh , Fh, 'replicate'));
	  RI_difcriBh  = abs(imfilter(  RI_criBh , Fh, 'replicate'));
	MLRI_difcriGrh = abs(imfilter(MLRI_criGrh, Fh, 'replicate'));
	MLRI_difcriGbh = abs(imfilter(MLRI_criGbh, Fh, 'replicate'));
	MLRI_difcriRh  = abs(imfilter(MLRI_criRh , Fh, 'replicate'));
	MLRI_difcriBh  = abs(imfilter(MLRI_criBh , Fh, 'replicate'));
	Fv = [-1;0;1];
	  RI_difcriGrv = abs(imfilter(  RI_criGrv, Fv, 'replicate'));
	  RI_difcriGbv = abs(imfilter(  RI_criGbv, Fv, 'replicate'));
	  RI_difcriRv  = abs(imfilter(  RI_criRv , Fv, 'replicate'));
	  RI_difcriBv  = abs(imfilter(  RI_criBv , Fv, 'replicate'));
	MLRI_difcriGrv = abs(imfilter(MLRI_criGrv, Fv, 'replicate'));
	MLRI_difcriGbv = abs(imfilter(MLRI_criGbv, Fv, 'replicate'));
	MLRI_difcriRv  = abs(imfilter(MLRI_criRv , Fv, 'replicate'));
	MLRI_difcriBv  = abs(imfilter(MLRI_criBv , Fv, 'replicate'));

	% absolute value of iteration criteria
	  RI_criGrh = abs(  RI_criGrh);
	  RI_criGbh = abs(  RI_criGbh);
	  RI_criRh  = abs(  RI_criRh );
	  RI_criBh  = abs(  RI_criBh );
	  RI_criGrv = abs(  RI_criGrv);
	  RI_criGbv = abs(  RI_criGbv);
	  RI_criRv  = abs(  RI_criRv );
	  RI_criBv  = abs(  RI_criBv );
	MLRI_criGrh = abs(MLRI_criGrh);
	MLRI_criGbh = abs(MLRI_criGbh);
	MLRI_criRh  = abs(MLRI_criRh );
	MLRI_criBh  = abs(MLRI_criBh );
	MLRI_criGrv = abs(MLRI_criGrv);
	MLRI_criGbv = abs(MLRI_criGbv);
	MLRI_criRv  = abs(MLRI_criRv );
	MLRI_criBv  = abs(MLRI_criBv );

    % add Gr and R (Gb and B) criteria residuals 
	  RI_criGRh = (  RI_criGrh +   RI_criRh).*Mrh;
	  RI_criGBh = (  RI_criGbh +   RI_criBh).*Mbh;
	  RI_criGRv = (  RI_criGrv +   RI_criRv).*Mrv;
	  RI_criGBv = (  RI_criGbv +   RI_criBv).*Mbv;
	MLRI_criGRh = (MLRI_criGrh + MLRI_criRh).*Mrh;
	MLRI_criGBh = (MLRI_criGbh + MLRI_criBh).*Mbh;
	MLRI_criGRv = (MLRI_criGrv + MLRI_criRv).*Mrv;
	MLRI_criGBv = (MLRI_criGbv + MLRI_criBv).*Mbv;

    % add Gr and R (Gb and B) gradient of criteria residuals	
	  RI_difcriGRh = (  RI_difcriGrh +   RI_difcriRh).*Mrh;
	  RI_difcriGBh = (  RI_difcriGbh +   RI_difcriBh).*Mbh;
	  RI_difcriGRv = (  RI_difcriGrv +   RI_difcriRv).*Mrv;
	  RI_difcriGBv = (  RI_difcriGbv +   RI_difcriBv).*Mbv;
	MLRI_difcriGRh = (MLRI_difcriGrh + MLRI_difcriRh).*Mrh;
	MLRI_difcriGBh = (MLRI_difcriGbh + MLRI_difcriBh).*Mbh;
	MLRI_difcriGRv = (MLRI_difcriGrv + MLRI_difcriRv).*Mrv;
	MLRI_difcriGBv = (MLRI_difcriGbv + MLRI_difcriBv).*Mbv;

	% directional map of iteration criteria
	  RI_crih =   RI_criGRh +   RI_criGBh;
	  RI_criv =   RI_criGRv +   RI_criGBv;
	MLRI_crih = MLRI_criGRh + MLRI_criGBh;
	MLRI_criv = MLRI_criGRv + MLRI_criGBv;
	
	% directional gradient map of iteration criteria
	  RI_difcrih =   RI_difcriGRh+  RI_difcriGBh;
	  RI_difcriv =   RI_difcriGRv+  RI_difcriGBv;
	MLRI_difcrih = MLRI_difcriGRh+MLRI_difcriGBh;
	MLRI_difcriv = MLRI_difcriGRv+MLRI_difcriGBv;

    % smoothing of iteration criteria
	sigma = 2;
	Fh = fspecial('gaussian', [5,5], sigma);
	  RI_crih = imfilter(  RI_crih, Fh, 'replicate');
	MLRI_crih = imfilter(MLRI_crih, Fh, 'replicate');
	  RI_difcrih = imfilter(  RI_difcrih, Fh, 'replicate');
	MLRI_difcrih = imfilter(MLRI_difcrih, Fh, 'replicate');
	Fv = fspecial('gaussian', [5,5],sigma);
	  RI_criv = imfilter(  RI_criv, Fv, 'replicate');
	MLRI_criv = imfilter(MLRI_criv, Fv, 'replicate');
	  RI_difcriv = imfilter(  RI_difcriv, Fv, 'replicate');
	MLRI_difcriv = imfilter(MLRI_difcriv, Fv, 'replicate');

    % calcualte iteration criteria (Eq.(13))
	  RI_wh = (  RI_crih.^2).*(  RI_difcrih);
	  RI_wv = (  RI_criv.^2).*(  RI_difcriv);
	MLRI_wh = (MLRI_crih.^2).*(MLRI_difcrih);
	MLRI_wv = (MLRI_criv.^2).*(MLRI_difcriv);

    % find smaller criteria pixels
	  RI_pih = find(   RI_wh<  RI_w2h );
	  RI_piv = find(   RI_wv<  RI_w2v );
	MLRI_pih = find( MLRI_wh<MLRI_w2h );
	MLRI_piv = find( MLRI_wv<MLRI_w2v );

    % guide updating
	  RI_Guidegrh = mosaic(:,:,2).*maskGr + RI_Grh; 
	  RI_Guidegbh = mosaic(:,:,2).*maskGb + RI_Gbh; 
	  RI_Guidegh  = RI_Guidegrh + RI_Guidegbh; 
	  RI_Guiderh  = mosaic(:,:,1) + RI_Rh;
	  RI_Guidebh  = mosaic(:,:,3) + RI_Bh;
	  RI_Guidegrv = mosaic(:,:,2).*maskGb + RI_Grv; 
	  RI_Guidegbv = mosaic(:,:,2).*maskGr + RI_Gbv; 
	  RI_Guidegv  = RI_Guidegrv + RI_Guidegbv; 
	  RI_Guiderv  = mosaic(:,:,1) + RI_Rv;
	  RI_Guidebv  = mosaic(:,:,3) + RI_Bv;
	MLRI_Guidegrh = mosaic(:,:,2).*maskGr + MLRI_Grh; 
	MLRI_Guidegbh = mosaic(:,:,2).*maskGb + MLRI_Gbh; 
	MLRI_Guidegh  = MLRI_Guidegrh + MLRI_Guidegbh; 
	MLRI_Guiderh  = mosaic(:,:,1) + MLRI_Rh;
	MLRI_Guidebh  = mosaic(:,:,3) + MLRI_Bh;
	MLRI_Guidegrv = mosaic(:,:,2).*maskGb + MLRI_Grv; 
	MLRI_Guidegbv = mosaic(:,:,2).*maskGr + MLRI_Gbv; 
	MLRI_Guidegv  = MLRI_Guidegrv + MLRI_Guidegbv; 
	MLRI_Guiderv  = mosaic(:,:,1) + MLRI_Rv;
	MLRI_Guidebv  = mosaic(:,:,3) + MLRI_Bv;

    % select smallest iteration criteria at each pixel (Eq.(14))
	  RI_Gh(  RI_pih) =   RI_Guidegh(  RI_pih);
	MLRI_Gh(MLRI_pih) = MLRI_Guidegh(MLRI_pih);
	  RI_Gv(  RI_piv) =   RI_Guidegv(  RI_piv);
	MLRI_Gv(MLRI_piv) = MLRI_Guidegv(MLRI_piv);

    % update minimum iteration criteria
	  RI_w2h(  RI_pih) =   RI_wh(  RI_pih);
	  RI_w2h(  RI_pih) =   RI_wh(  RI_pih);
	  RI_w2v(  RI_piv) =   RI_wv(  RI_piv);
	  RI_w2v(  RI_piv) =   RI_wv(  RI_piv);
	MLRI_w2h(MLRI_pih) = MLRI_wh(MLRI_pih);
	MLRI_w2h(MLRI_pih) = MLRI_wh(MLRI_pih);
	MLRI_w2v(MLRI_piv) = MLRI_wv(MLRI_piv);
	MLRI_w2v(MLRI_piv) = MLRI_wv(MLRI_piv);

	% guided filter window size update
	h = h+1;
	v = v+1;
	h2 = h2+1;
	v2 = v2+1;
end

%% Step(iii): adaptive combining
% combining weight (Eq.(16))
  RI_w2h = 1./(  RI_w2h + 1e-10);
  RI_w2v = 1./(  RI_w2v + 1e-10);
MLRI_w2h = 1./(MLRI_w2h + 1e-10);
MLRI_w2v = 1./(MLRI_w2v + 1e-10);
w = RI_w2h + RI_w2v + MLRI_w2h + MLRI_w2v;

% combining (Eq.(15))
green = (RI_w2h.*RI_Gh + RI_w2v.*RI_Gv + MLRI_w2h.*MLRI_Gh + MLRI_w2v.*MLRI_Gv)./(w+1e-32);

% final output
green = green.*imask(:,:,2) + mosaic(:,:,2);
green = clip(green, 0, DynamicRange);

end

