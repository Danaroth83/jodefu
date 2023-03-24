function [ methods_list ] = load_methodslist( methodsclass,flag_removeEXP )
%LOAD_METHODSLIST Loads a list of fusion methods to perform
%   Loads a particular list of fusion methods to perform, according to a
%   specific label; it also allows to input a specific method or list of
%   methods
%   if flag_removeEXP=0, it removes the straight interpolation method

methods_list={'EXP','PCA','IHS','Brovey','BDSD','GS','GSA','GSA_local','PRACS','HPF','SFIM',...
    'Indusion','ATWT','AWLP','ATWT_M2','ATWT_M3','ATWT_mod','MTF_GLP','MTF_GLP_HPM_PP',...
    'MF_HG','MTF_GLP_HPM','MTF_GLP_CBD','GFPCA','CNMF','BayesNaive_FE','BayesSparse_FE',...
    'HySure','G_BDSD','L_BDSD','C_BDSD','BayesNaive','BayesSparse','BayesML',...
    'BayesNLM','none'};
%,'DMSC','DMSC_period','DMSC_Bayes','DMSC_mindis','DMSC_random',...
%    'DMSC_weighted','DMSC_weighted2','DMSC_weighted3','DMSC_weighted4',...
%    'CASSI','CASSI_TV','CASSI_TV2','CFA_WAV'};

if nargin<2, flag_removeEXP=0; end


if iscell(methodsclass)
    idx_methods=[];
    for ii2=1:numel(methodsclass)
        for ii=1:numel(methods_list)
            if strcmpi(strrep(methodsclass{ii2},'-','_'),methods_list{ii})
                idx_methods=cat(1,idx_methods,ii);
            end
        end
    end
elseif strncmpi(methodsclass,'all',3)
    idx_methods=1:numel(methods_list);
elseif strncmpi(methodsclass,'IGARSS2016',10)
    idx_methods=[1:2,7,14,21,24,27];
elseif strncmpi(methodsclass,'BioucasDias',11)
    idx_methods=[1,6:7,11,17,21,23:27];
elseif strncmpi(methodsclass,'MTFGLP',6)
    idx_methods=[1,21:22];
elseif strncmpi(methodsclass,'removepersonal',14)
    idx_methods=[1:7,9:15,18:27];
elseif strncmpi(methodsclass,'basic',5)
    idx_methods=[1:4,6:7,13,21:22];
elseif strncmpi(methodsclass,'removeslow',10)
    idx_methods=[1:25,28:30];
elseif strncmpi(methodsclass,'originalToolbox',15) || strncmpi(methodsclass,'Alparone',8)
    idx_methods=[1:7,9:22];
elseif strncmpi(methodsclass,'usualtest',9)
    idx_methods=[1:4,6:7,11:14,17,21:22,24:25,27:28];
elseif strncmpi(methodsclass,'WHISPERS2016',12)
    idx_methods=[1,11,13,20,18,21,22];
% elseif strncmpi(methodsclass,'custom',6)
%     setup_flags;
%     idx_methods=1:numel(methods_list);
%     for ii=1:numel(methods_list)
%         if eval(['flag_',methods_list{ii}])==0
%             idx_methods=idx_methods(idx_methods~=ii);
%         end
%     end
else
    for ii=1:numel(methods_list)
        if strcmpi(strrep(methodsclass,'-','_'),methods_list{ii})
            idx_methods=ii;
        end
    end
end  

if flag_removeEXP==1, idx_methods=idx_methods(idx_methods~=1); end
methods_list=methods_list(idx_methods);

end

