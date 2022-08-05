function [gx,gy] = sgrad(u,n)
  %%----------------------------- SGRAD -----------------------------
  %
  %                        Shannon gradient
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
  % Usage: [gx,gy] = sgrad(u,n)
  %
  % + u : input image (2D matrix)
  % + n : oversampling parameter (use STVn regularizer)
  % + gx : Shannon gradient along X-axis (horizontal)
  % + gy : Shannon gradient along Y-axis (vertical)
  %
  % v1.0 (12/2016): initial version.
  
  %% Initialization
  [N,M] = size(u); % size of the image domain Omega (M = width, N=height).
  nM = n*M; nN = n*N; % width & height of the refined domain Omega_n.
  M_is_even = (mod(M,2)==0);
  N_is_even = (mod(N,2)==0);

  %% dispatch dft_u entries over a symmetrized grid 
  a = -floor(M/2)+(0:M-1);
  b = -floor(N/2)+(0:N-1);
  if(M_is_even); a = [a,M/2]; end
  if(N_is_even); b = [b,N/2]; end
  a = ifftshift(a);
  b = ifftshift(b);
  [A,B] = meshgrid(a,b);
  cof = ones(size(A));
  if (M_is_even); cof(abs(A) == M/2) = 0.5; end
  if (N_is_even); cof(abs(B) == N/2) = 0.5; end
  if (M_is_even && N_is_even); cof((abs(A) == M/2)&(abs(B) == N/2)) = 0.25; end
  dft_u = fft2(u);
  dft_u_augm = cof .* dft_u(1+mod(N+b,N),1+mod(M+a,M));

  %% compute (dft_gx,dft_gy) = DFT(Shannon gradient of u) 
  dft_gx = zeros(nN,nM);
  dft_gy = zeros(nN,nM);
  Omega_a = 1+mod(nM+a,nM); % location of frequencies [-M/2,M/2] in dft_gx and dft_gy along the horizontal axis
  Omega_b = 1+mod(nN+b,nN); % location of frequencies [-N/2,N/2] in dft_gx and dft_gy along the vertical axis
  dft_gx(Omega_b,Omega_a) = n^2*2*1i*pi * (dft_u_augm.*A)/M;
  dft_gy(Omega_b,Omega_a) = n^2*2*1i*pi * (dft_u_augm.*B)/N;

  %% compute inverse DFT to obtain (gx,gy) = Shannon gradient of u
  gx = real(ifft2(dft_gx)); % nonzero imaginary part is only due to numerical errors (max(abs(imag(gx))) is close to the machine epsilon)
  gy = real(ifft2(dft_gy)); % nonzero imaginary part is only due to numerical errors (max(abs(imag(gy))) is close to the machine epsilon)

end

