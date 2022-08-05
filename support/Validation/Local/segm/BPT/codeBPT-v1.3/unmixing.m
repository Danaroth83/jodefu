function model = unmixing(pixels,method)

[nbands,nsamples] = size(pixels);
model = struct;
% Induce the endmembers
runs = 20;
[w,Rn] = estNoise(pixels,'additive');
[p_hys,~]=hysime(pixels,w,Rn);
if p_hys > 0 && p_hys < nsamples
    results = EIA_1D(pixels,'vca','p',p_hys,'runs',runs);
    % Obtain the optimum E in base to RMSE
    RMSE_opt = [];
    it_opt = 1;
    A_opt = [];
    RMSE_avg_opt = [];
    RMSE_max_opt = [];
    for i=1:runs
        fprintf('Iteration %g\n',i);
        % Endmembers
        E = results(i,1).E;
        % FCLSU
        A = FCLSU(pixels,E);
        % Reconstruction
        H = E*A';
        % RMSE
        RMSE = sqrt((sum((pixels - H).^2))./nbands);
        RMSE_avg = mean(RMSE);
        RMSE_max = max(RMSE);
        if strcmp(method,'max')
            if RMSE_max < RMSE_max_opt
                RMSE_opt = RMSE;
                RMSE_avg_opt = RMSE_avg;
                RMSE_max_opt = RMSE_max;
                it_opt = i;
                A_opt = A;
            end
        else
            if RMSE_avg < RMSE_avg_opt
                RMSE_opt = RMSE;
                RMSE_avg_opt = RMSE_avg;
                RMSE_max_opt = RMSE_max;
                it_opt = i;
                A_opt = A;
            end
        end
    end
    E_opt = results(it_opt,1).E;
    p_opt = size(E,2);
    C_opt = results(it_opt,1).C;
    model.method = 'vca';
else
    E_opt = mean(pixels,2);
    C_opt = [];
    p_opt = 1;
    % FCLSU
    A_opt = FCLSU(pixels,E_opt);
    % Reconstruction
    H = E_opt*A_opt';
    % SSE
    RMSE_opt = (sqrt(sum((pixels - H).^2))./nbands);
    RMSE_avg_opt = mean(RMSE_opt);
    RMSE_max_opt = max(RMSE_opt);
    model.method = 'mean';
end
model.E = E_opt;
model.C = C_opt;
model.p = p_opt;
model.A = A_opt;
model.RMSE = RMSE_opt;
model.RMSE_avg = RMSE_avg_opt;
model.RMSE_max = RMSE_max_opt;