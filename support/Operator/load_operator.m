%% Load Operator
%
% Author:
% Daniele Picone
%
% Description:
% Function wrapper to load various kind of operators, generally intended to
% be applied to color images.
%
% Usage:
% Op=load_operator(choice);
%
% Input:
% choice: String defining the desired operator, among the following
%         categories
%       'WAVn_???': 2D Discrete Wavelet transform, where 'n' is the amount 
%                   of decomposition levels and '???' is a string defining 
%                   the wavelet type (some possible choices: 'db4','sym8',
%                   'haar','coif2' and all matLAB supported)  Eg: 'WAV3_db8'
%       'CASn_???': Same as above but Discrete Cosine Transform is also
%                   performed over the third dimension
%       'TV_?'    : Total Variation operator, where '?' is a char
%                   describing the type
%              - 'TV_c': Classic Total Variation
%              - 'TV_u': Upwind Total Variation [1]
%              - 'TV_p': New discrete Total Variation [2]
%              - 'TV_m': Variant of the above
%              - 'TV_h': Hessian
%              - 'TV_g': Total Generalized Variation (TGV) [3]
%              - 'TV_s': Shannon Total Variation [4]
%       'NORM_???': Norm operator for vector quantities [5,6] (possible 
%                   choices: 'l1','l2','l211','l221','linf11', 'linfinf1,
%                   'linf21','l2inf1','S1l1', 'S1linf')
%        'none'   : Returns the input as is
% sizes: The sizes of the input image (just necessary for wavelet transforms)
%
% Output:
% Op: Struct defining some properties of the desired operator f(x); its fields
%     depend on whchi option is available for each operator and include:
%       direct  : The function handle for the direct operation f(x)
%       adjoint : The function handle of the adjoint of the above 
%                 (just for linear operators)
%       inverse : The function handle for the inverse operation f^(-1)(x)
%       norm    : The operator norm 
%       prox    : The function handle for the proximal operator (just for functions in R+)
%                 It has two inputs (x,gamma) for the operator of gamma f(x)
%       proxconj: The function handle for the proximal operator of the 
%                 (Fenchel) convex conjugate f*(x); it has in general
%                 three inputs (x,gamma,lambda) for the function 
%                 gamma [lambda f(x)]* (just available for functions in R+)
%
% References:
% [1] Chambolle A., Levine S. E. and Lucier B. J., "An Upwind 
% Finite-Difference Method for Total Variation–Based Image Smoothing" SIAM 
% Journal on Imaging Sciences, 2011, 4, 277-299
% [2] Condat L., "Discrete Total Variation: New Definition and 
% Minimization" SIAM Journal on Imaging Sciences, 2017, 10, 1258-1290
% [3] Bredies K., Kunisch K. and Pock T. "Total generalized variation",
% SIAM Journal on Imaging Sciences, 2010, 3, 492-526
% [4] Abergel R. and Moisan L., "The Shannon Total Variation", Journal of 
% Mathematical Imaging and Vision, Springer Nature, 2017, 59, 341-370
% [5] Duran J., Moeller M., Sbert C. and Cremers D., "Collaborative total 
% variation: a general framework for vectorial TV models", SIAM Journal on 
% Imaging Sciences, 2016, 9, 116-151
% [6] Duran J., Moeller M., Sbert C. and Cremers D., "On the implementation
% of collaborative TV regularization: Application to cartoon+ texture
% decomposition" IPOL, 2016, 6, 27-74

function Op=load_operator(choice,sizes)

    current_folder=fileparts(mfilename('fullpath'));

    if strncmpi(choice,'wav',3) || strncmpi(choice,'cas',3)
        
        addpath(fullfile(current_folder,'Wavelet'));
        
        sep=strfind(choice,'_');
        wavlevels=str2double(choice(4:sep-1));
        wav_type=choice(sep+1:end);

        % % If the wavelet toolbox is available
        % wav_type=choice(sep+1:end);
        % [LoD,HiD,LoR,HiR] = wfilters(wav_type);

        % Bypass for who doesn't own Matlab Wavelet toolbox
        addpath(fullfile(current_folder,'..','Matlab_toolboxes','wavelet'));
        if strcmpi(wav_type,'sym8')
            W = [0.00133639669640  -0.00021419715012  -0.01057284326418   0.00269319437688  ...
                 0.03474523295559  -0.01924676063167  -0.03673125438038   0.25769933518654  ...
                 0.54955331526901   0.34037267359439  -0.04332680770282  -0.10132432764282  ...
                 0.00537930587524   0.02241181152181  -0.00038334544811  -0.00239172925575];
        elseif any(strcmpi(wav_type,{'haar','db1'}))
            W = [0.5,0.5];
        elseif strcmpi(wav_type,'db4')
            W = [ 0.16290171402562   0.50547285754565   0.44610006912319  -0.01978751311791 ...
                 -0.13225358368437   0.02180815023739   0.02325180053556  -0.00749349466513];
        elseif strcmpi(wav_type,'db8')
            W = [ 0.03847781105406   0.22123362357624   0.47774307521438   0.41390826621166 ...
                 -0.01119286766665  -0.20082931639111   0.00033409704628   0.09103817842345 ...
                 -0.01228195052300  -0.03117510332533   0.00988607964808   0.00618442240954 ...
                 -0.00344385962813  -0.00027700227421   0.00047761485533  -0.00008306863060];
        elseif strcmpi(wav_type,'coif1')
            W = [-0.05142972847100   0.23892972847100   0.60285945694200   0.27214054305800 ...
                 -0.05142972847100  -0.01107027152900];
        elseif strcmpi(wav_type,'coif2')
            W = [ 0.01158759673900  -0.02932013798000  -0.04763959031000   0.27302104653500 ...
                  0.57468239385700   0.29486719369600  -0.05408560709200  -0.04202648046100 ...
                  0.01674441016300   0.00396788361300  -0.00128920335600  -0.00050950539900];
        else
            error('Just Symlet 8 are currently implemented');
        end
        W = W/sum(W); LoR = sqrt(2)*W; HiR = LoR(end:-1:1);
        p=0; first = 2-rem(p,2);
        HiR(first:2:end) = -HiR(first:2:end); HiD = flip(HiR); LoD = flip(LoR);

        W1 = @(x) ipermute(idct(permute(x,[3,1,2])),[3,1,2]);
        W2 = @(x) idwt2_custom(x,LoR,HiR,wavlevels,sizes);
        W1I = @(x) ipermute(dct(permute(x,[3,1,2])),[3,1,2]);
        W2I = @(x) dwt2_custom(x,LoD,HiD,wavlevels);
        
        if strcmpi(choice(1),'c')
            opdirect  = @(x) W1(W2(x));
            opinverse = @(x) W2I(W1I(x));
        else
            opdirect  = @(x) W2(x);
            opinverse = @(x) W2I(x);
        end
        opadjoint=opinverse;

        Op.direct  = @(x) reshape(opdirect(reshape(x,size(x,1),size(x,2),size(x,3),[])),size(x));
        Op.adjoint = @(x) reshape(opinverse(reshape(x,size(x,1),size(x,2),size(x,3),[])),size(x));
        Op.inverse = @(x) reshape(opadjoint(reshape(x,size(x,1),size(x,2),size(x,3),[])),size(x));
        Op.norm    = 1; 

    elseif strncmpi(choice,'TV',2)

        addpath(fullfile(current_folder,'Total_Variation'));
        sep=strfind(choice,'_');
        Op=TV_operator(choice(sep+1:end));
        
    elseif strncmpi(choice,'Norm',4)
        
        addpath(fullfile(current_folder,'Norm'));
        sep=strfind(choice,'_');
        Op=CTV_operator('norm',choice(sep+1:end));
    
    elseif strcmpi(choice,'none')
        
        Op.direct   = @(x) x;
        Op.adjoint  = @(x) x;
        Op.inverse  = @(x) x;
        Op.prox     = @(x,gamma) x;
        Op.proxconj = @(x,gamma,lambda) zeros(size(x));
        Op.norm     = 1;
        
    elseif isnumeric(choice)
        
        % IT supposes it is a matrix multiplication
        
        A=choice;
            
        Op.direct   = @(x) A*x;
        Op.adjoint  = @(u) A'*u;
        Op.inverse  = @(x) A\x;
        Op.norm     = max(svd(A));
    
    else
        error('Choice is unknown');     
    end
end
