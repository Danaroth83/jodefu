function div = sdiv(px,py,n)

%%                      Shannon divergence
%
% Authors: Rémy Abergel & Lionel Moisan.
%
%   This program is freely available on the web page
%
%   http://www.math-info.univ-paris5.fr/~rabergel/
%
%   We hope that you will find it useful. If you use it for a
%   publication, please mention this web page and the paper it makes
%   reference to. If you modify this program, please indicate in the
%   code that you did so and leave this message. You can also report
%   bugs or suggestions to Remy.Abergel [AT] gmail.com and
%   Lionel.Moisan [AT] parisdescartes.fr
%
% Usage: div = sdiv(px,py,n)
%
% + n : oversampling parameter (use STVn regularizer)
% + px : first component of the vector field p = (px,py)
% + py : second component of the vector field p = (px,py)
%
% v1.0 (12/2016): initial Scilab version
%

%% Initialization 
[nN,nM] = size(px);
M = nM/n; N = nN/n;
M_is_even = (mod(M,2)==0);
N_is_even = (mod(N,2)==0);
dft_div = zeros(N,M);
dft_px = fft2(px);
dft_py = fft2(py);

%% compute frequency grids 
a0 = -floor((M-1)/2) : floor((M-1)/2);
b0 = -floor((N-1)/2) : floor((N-1)/2);
[A0,B0] = meshgrid(a0,b0);
Omega_a0 = 1+mod(M+a0,M); % location of frequencies [-M/2,M/2] in dft_div along the horizontal axis
Omega_b0 = 1+mod(N+b0,N); % location of frequencies [-N/2,N/2] in dft_div along the vertical axis
Omega_aa0 = 1+mod(nM+a0,nM); % location of frequencies [-M/2,M/2] in dft_px & dft_py along the horizontal axis
Omega_bb0 = 1+mod(nN+b0,nN); % location of frequencies [-N/2,N/2] in dft_px & dft_py along the vertical axis

%% compute div = Shannon divergence of p:=(px,py) 
dft_div(Omega_b0,Omega_a0) = 2*1i*pi*(dft_px(Omega_bb0,Omega_aa0).*(A0/M) + dft_py(Omega_bb0,Omega_aa0).*(B0/N));
if (M_is_even)
    a = -M/2;
    adr_a_plus = 1+mod(nM+a,nM); % location of frequency a in dft_px & dft_py along the horizontal axis
    adr_a_minus = 1+mod(nM-a,nM); % location of frequency -a in dft_px & dft_py along the horizontal axis
    hpx = 0.5*(dft_px(Omega_bb0,adr_a_plus) - dft_px(Omega_bb0,adr_a_minus));
    hpy = 0.5*(dft_py(Omega_bb0,adr_a_plus) + dft_py(Omega_bb0,adr_a_minus));
    dft_div(Omega_b0,1+M/2) = 2*1i*pi*((a/M)*hpx + (b0'/N).*hpy);
end
if (N_is_even) 
    b = -N/2;
    adr_b_plus = 1+mod(nN+b,nN); % location of frequency b in dft_px & dft_py along the vertical axis
    adr_b_minus = 1+mod(nN-b,nN); % location of frequency -b in dft_px & dft_py along the vertical axis
    hpx = 0.5*(dft_px(adr_b_plus,Omega_aa0) + dft_px(adr_b_minus,Omega_aa0));
    hpy = 0.5*(dft_py(adr_b_plus,Omega_aa0) - dft_py(adr_b_minus,Omega_aa0));
    dft_div(1+N/2,Omega_a0) = 2*1i*pi*((a0/M).*hpx + (b/N)*hpy);
end
if (M_is_even) && (N_is_even) 
    a = -M/2;
    b = -N/2;
    adr_a_plus = 1+mod(nM+a,nM); % location of frequency a in dft_px & dft_py along the horizontal axis
    adr_a_minus = 1+mod(nM-a,nM); % location of frequency -a in dft_px & dft_py along the horizontal axis
    adr_b_plus = 1+mod(nN+b,nN); % location of frequency b in dft_px & dft_py along the vertical axis
    adr_b_minus = 1+mod(nN-b,nN); % location of frequency -b in dft_px & dft_py along the vertical axis
    hpx = 0.25*(dft_px(adr_b_plus,adr_a_plus) - dft_px(adr_b_plus,adr_a_minus) + dft_px(adr_b_minus,adr_a_plus) - dft_px(adr_b_minus,adr_a_minus));
    hpy = 0.25*(dft_py(adr_b_plus,adr_a_plus) + dft_py(adr_b_plus,adr_a_minus) - dft_py(adr_b_minus,adr_a_plus) - dft_py(adr_b_minus,adr_a_minus));
    dft_div(1+N/2,1+M/2) = 2*1i*pi*((a/M)*hpx + (b/N)*hpy);
end

%% compute inverse DFT to obtain (gx,gy) = Shannon gradient of u 
div = real(ifft2(dft_div)); % nonzero imaginary part is only due to numerical errors (max(abs(imag(div))) is close to the machine epsilon)

end

