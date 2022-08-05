function [ I_PAN_LR, I_PAN_eq ] = PAN_equalization( I_PAN, I_MS, flag_degradation, flag_equalization, ratio, GNyq, GNyq_PAN, I_GT )
%PAN_EQUALIZATION Summary of this function goes here
%   Detailed explanation goes here

current_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'..','scale'));

if nargin<=2, flag_degradation=0; end
if nargin<=3, flag_equalization=flag_degradation; end
if nargin<=4 && flag_degradation~=0, error('scale ratio was not provided'); end
if nargin<=7 && flag_equalization==1, error('GT was not provided'); end

if flag_degradation==0
    I_PAN_LR=I_PAN;
elseif flag_degradation==1
    I_PAN_LR=imresize(imresize(I_PAN,1/ratio),ratio);
elseif flag_degradation==2
    I_PAN_LR=resize_images_PAN(I_PAN,ratio,GNyq_PAN,0,0);
elseif flag_degradation==3
    I_PAN_LR=LowPass_Bspline(I_PAN,ratio,3);
elseif flag_degradation==4
    fcut=1/ratio;
    N = 21;
    [f1,f2] = freqspace(N,'meshgrid');
    Hd = ones(N);
    Hd((abs(f1)>fcut)|(abs(f2)>fcut)) = 0;
    h = fwind1(Hd,hamming(N));
    I_PAN_LR = imfilter(I_PAN,real(h),'replicate');
elseif flag_degradation==5
    I_PAN_LR=resize_images_MS(repmat(I_PAN,[1,1,size(I_MS,3)]),ratio,GNyq,0,0);
    % I_PAN_LR=sum(I_PAN_LR,3)/size(I_PAN_LR,3);
end

if nargout==2
    if flag_equalization==0 || flag_equalization==1
        I_PAN_load=I_PAN_LR;
        if flag_equalization==0
            I_MS_load=I_MS;
        elseif flag_equalization==1
            I_MS_load=I_GT;
        end

        [L1,L2,Lb]=size(I_PAN_load);
        Dim=L1*L2;
        I_PAN_temp=reshape(I_PAN_load,[Dim,Lb]);
        mean2_PAN=sum(I_PAN_temp,1)/Dim;
        std2_PAN=sqrt(sum(I_PAN_temp.^2,1)/(Dim-1)-mean2_PAN.^2);

        [N1,N2,Nb]=size(I_MS_load);
        Dim2=N1*N2;
        I_MS_temp=reshape(I_MS_load,[Dim2,Nb]);
        mean2_MS=sum(I_MS_temp,1)/Dim2;
        std2_MS=sqrt(sum(I_MS_temp.^2,1)/(Dim2-1)-mean2_MS.^2);

        I_PAN_eq=(repmat(I_PAN,[1,1,Lb])-repmat(reshape(mean2_PAN,[1,1,Lb]),[L1,L2,1]))./repmat(reshape(std2_PAN,[1,1,Lb]),[L1,L2,1]);
        I_PAN_eq=repmat(I_PAN_eq,[1,1,Nb-Lb+1]).*repmat(reshape(std2_MS,[1,1,Nb]),[L1,L2,1])+repmat(reshape(mean2_MS,[1,1,Nb]),[L1,L2,1]);
    elseif flag_equalization==2
        [L1,L2,Lb]=size(I_PAN);
        Dim=L1*L2;
        I_PAN_temp=reshape(I_PAN,[Dim,Lb]);
        mean2_PAN=sum(I_PAN_temp,1)/Dim;
        % std2_PAN=sqrt(sum(I_PAN_temp.^2,1)/(Dim-1)-mean2_PAN.^2);
        
        [N1,N2,Nb]=size(I_MS);
        Dim2=N1*N2;
        I_MS_v=reshape(I_MS,[Dim2,Nb]);
        mean2_MS=sum(I_MS_v,1)/Dim2;
        std2_MS=sqrt(sum(I_MS_v.^2,1)/(Dim2-1)-mean2_MS.^2);
        alpha_MS=estimation_alpha2(I_MS,sum(I_PAN_LR,3)/size(I_PAN_LR,3),'global');
        std2_den=sqrt(sum((alpha_MS(1:end-1)'.^2).*(std2_MS.^2)));
        
        I_PAN_eq=(I_PAN-repmat(reshape(mean2_PAN,[1,1,Lb]),[L1,L2,1]))./repmat(reshape(std2_den,[1,1,Lb]),[L1,L2,1]);
        I_PAN_eq=repmat(I_PAN_eq,[1,1,Nb-Lb+1]).*repmat(reshape(std2_MS,[1,1,Nb]),[L1,L2,1])+repmat(reshape(mean2_MS,[1,1,Nb]),[L1,L2,1]);
    elseif flag_equalization==3
        [L1,L2,Lb]=size(I_PAN);
        Dim=L1*L2;
        I_PAN_temp=reshape(I_PAN,[Dim,Lb]);
        mean2_PAN=sum(I_PAN_temp,1)/Dim;
        % std2_PAN=sqrt(sum(I_PAN_temp.^2,1)/(Dim-1)-mean2_PAN.^2);
        
        [N1,N2,Nb]=size(I_MS);
        Dim2=N1*N2;
        I_MS_v=reshape(I_MS,[Dim2,Nb]);
        [~,I_MS_temp]=pca(I_MS_v);
        I_MS_PCA=reshape(I_MS_temp,[N1,N2,Nb]);
        alpha_PCA=estimation_alpha(I_MS_PCA,sum(I_PAN_LR,3)/size(I_PAN_LR,3),'global');
        mean2_PCA=sum(I_MS_temp,1)/Dim2;
        std2_PCA=sqrt(sum(I_MS_temp.^2,1)/(Dim2-1)-mean2_PCA.^2);
        std2_den=sqrt(sum((alpha_PCA'.^2).*(std2_PCA.^2)));
        mean2_MS=sum(I_MS_v,1)/Dim2;
        std2_MS=sqrt(sum(I_MS_v.^2,1)/(Dim2-1)-mean2_MS.^2);
        
        I_PAN_eq=(I_PAN-repmat(reshape(mean2_PAN,[1,1,Lb]),[L1,L2,1]))./repmat(reshape(std2_den,[1,1,Lb]),[L1,L2,1]);
        I_PAN_eq=repmat(I_PAN_eq,[1,1,Nb-Lb+1]).*repmat(reshape(std2_MS,[1,1,Nb]),[L1,L2,1])+repmat(reshape(mean2_MS,[1,1,Nb]),[L1,L2,1]);
    end
        
end
end