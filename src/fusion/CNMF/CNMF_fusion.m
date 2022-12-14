function Out = CNMF_fusion(HSI,MSI,mask)
%--------------------------------------------------------------------------
% COUPLED NONNEGATIVE MATRIX FACTORIZATION (CNMF)
%
% USAGE
%       Out = CNMF_fusion(HSI,MSI,mask)
%
% INPUT
%       HSI : Low-spatial-resolution HS image (rows2,cols2,bands2)
%       MSI : MS image (rows1,cols1,bands1)
%       mask: (optional) Binary mask for processing (rows2,cols2)
%
% OUTPUT
%       Out   : High-spatial-resolution HS image (rows1,cols1,bands2)
%
% Copyright (c) 2015 Naoto Yokoya All rights reserved.
% Email: yokoya@sal.rcast.u-tokyo.ac.jp
% Update: 2015/02/10
%
% REFERENCES
% [1] N. Yokoya, T. Yairi, and A. Iwasaki, "Coupled nonnegative matrix 
% factorization unmixing for hyperspectral and multispectral data fusion," 
% IEEE Trans. Geosci. Remote Sens., vol. 50, no. 2, pp. 528-537, 2012.
% [2] N. Yokoya, T. Yairi, and A. Iwasaki, "Hyperspectral, multispectral, 
% and panchromatic data fusion based on non-negative matrix factorization," 
% Proc. WHISPERS, Lisbon, Portugal, Jun. 6-9, 2011.
% [3] N. Yokoya, N. Mayumi, and A. Iwasaki, "Cross-calibration for data 
% fusion of EO-1/Hyperion and Terra/ASTER," IEEE J. Sel. Topics Appl. 
% Earth Observ. Remote Sens., vol. 6, no. 2, pp. 419-426, 2013.
%--------------------------------------------------------------------------

verbose = 'off'; % default

% masking mode
if nargin == 2
    masking = 0;
elseif nargin == 3
    masking = 1;
    mask2 = imresize(mask,2,'nearest');
else
    disp('Please check the usage of CNMF_fusion.m');
end

% image size
[rows1,cols1,bands1] = size(MSI);
[rows2,cols2,bands2] = size(HSI);

% estimation of R
if masking == 1
    R = estR(HSI,MSI,mask);
    for b = 1:bands1
        msi = reshape(MSI(:,:,b),rows1,cols1);
        msi = msi - R(b,end); % modify bias of msi
        msi(msi<0) = 0;
        MSI(:,:,b) = msi;
    end
    R = R(:,1:end-1);
else
    R = estR(HSI,MSI);
end

% parameters
th_h = 1e-8; % Threshold of change ratio in inner loop for HS unmixing
th_m = 1e-8; % Threshold of change ratio in inner loop for MS unmixing
th2 = 1e-2; % Threshold of change ratio in outer loop
sum2one = 2*(mean(MSI(:))/0.7455)^0.5/bands1^3; % Parameter of sum-to-one constraint
if bands1 == 1
    I1 = 75; % Maximum iteration of inner loop
    I2 = 1; % Maximum iteration of outer loop
else
    % For synthetic data: I1 = 100-300, I2 = 3-5
    % For real data     : I1 = 100-300, I2 = 1-2
    I1 = 100; % Maximum iteration of inner loop
    I2 = 1; % Maximum iteration of outer loop
end

% avoid nonnegative values
HSI(HSI<0) = 0;
MSI(MSI<0) = 0;

w = rows1/rows2;

switch masking
    case 0
        HSI = reshape(HSI,[],bands2)';
        MSI = reshape(MSI,[],bands1)';
    case 1
        HSI = reshape(HSI,[],bands2);
        MSI = reshape(MSI,[],bands1);
        % masking
        HSI = HSI(mask,:)';
        MSI = MSI(mask2,:)';
end

% number of endmembers
M = vd(HSI,5*10^-2); % M can be automatically defined, for example, by VD
%M = 20; % predefined
if strcmp (verbose, 'on'), disp(['number of endmembers: ' num2str(M)]); end

%% Initializatioin
switch masking
    case 0
        Out = CNMF_init(rows1,cols1,w,M,HSI,MSI,sum2one,I1,th_h,th_m,R);
    case 1
        Out = CNMF_init(rows1,cols1,w,M,HSI,MSI,sum2one,I1,th_h,th_m,R,mask);
end
HSI = Out.hyper;
MSI = Out.multi;
W_hyper = Out.W_hyper;
H_hyper = Out.H_hyper;
W_multi = Out.W_multi;
H_multi = Out.H_multi;
cost(1,1) = Out.RMSE_h;
cost(2,1) = Out.RMSE_m;

%% Iteration
for i = 1:I2
    switch masking
        case 0
            Out = CNMF_ite(rows1,cols1,w,M,HSI,MSI,W_hyper,H_hyper,W_multi,H_multi,I1,th_h,th_m,I2,i,R);
        case 1
            Out = CNMF_ite(rows1,cols1,w,M,HSI,MSI,W_hyper,H_hyper,W_multi,H_multi,I1,th_h,th_m,I2,i,R,mask);
    end
    cost(1,i+1) = Out.RMSE_h;
    if I2 > 1
        cost(2,i+1) = Out.RMSE_m;
    end
    if (cost(1,i)-cost(1,i+1))/cost(1,i)>th2 && (cost(2,i)-cost(2,i+1))/cost(2,i)>th2 && i<I2
        W_hyper=Out.W_hyper;
        H_hyper=Out.H_hyper;
        W_multi=Out.W_multi2;
        H_multi=Out.H_multi2;
    elseif i==I2
        W_hyper=Out.W_hyper;
        H_hyper=Out.H_hyper;
    else
        H_multi=Out.H_multi;
        break;
    end
end

switch masking
    case 0
        Out = reshape((W_hyper(1:bands2,:)*H_multi)',rows1,cols1,bands2);
    case 1
        Out = zeros(rows1*cols1,bands2);
        Out(mask2,:) = (W_hyper(1:bands2,:)*H_multi)';
        Out = reshape(Out,rows1,cols1,bands2);
end
