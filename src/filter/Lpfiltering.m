function h_out=Lpfiltering(GNyq,ratio,N)

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'Matlab_toolboxes'));

if nargin<=2, N = 41; end
if length(ratio)==1, ratio=[ratio,ratio]; end
fcut = 1 ./ ratio;
Nb = length(GNyq);

h_out=[];

for ii = 1 : Nb
    alpha = sqrt(((N-1)*(fcut/2)).^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', [N,1], alpha(2))*fspecial('gaussian', [1,N], alpha(1));
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    h = real(h);
    % I_MS_LR(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'replicate');
    h_out=cat(3,h_out,h);
end