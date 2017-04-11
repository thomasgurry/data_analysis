% Script to extract all outputs from cluster submission.
% Input total number of simulations, and number of parameters.

nsimulations = 1296;
nparams = 4;

str1 = 'output';
str2 = '.mat';

ncompleted = 0;

paramVals = zeros(nsimulations,nparams);
fvals = zeros(nsimulations,1);

for i = 1:nsimulations
    
    fullstr = strcat(str1,num2str(i),str2);
    
    fid = fopen(fullstr);
    
    if(fid ~= -1)
        
        ncompleted = ncompleted + 1;
        tmp = load(fullstr);
        paramVals(ncompleted,:) = tmp.out{1};
        fvals(ncompleted) = tmp.out{2};
        fclose(fid);
        
    end
    
end
        
paramVals = paramVals(1:ncompleted,:);
fvals = fvals(1:ncompleted);

% Sort the results by objective function value

indices = zeros(ncompleted,1);
sortedFvals = sort(fvals);

for i = 1:ncompleted
    
    indices(i) = min(find(fvals == sortedFvals(i)));
    
end

sortedParams = paramVals(indices,:);



%% Calculate observables for parameters resulting from optimisation. 
%% ONLY APPLICABLE TO NEGATIVE FEEDBACK.

all_observables = {};

% labels

all_observables{1} = {'RNA level','Protein level','RNA variance','Protein variance','Covariance','RNA CV','Protein CV'};

%% First: parameters effect on SGE

observables = zeros(4,7);

paramsSGE = [10,10,1,0.001,200];
paramsNF1000 = [10,0.001,0.0026,0.001,200];
paramsNF200 = [10,10,0.4,0.001,200];
paramsNF50 = [10,10,0.0531,0.001,200];

[x] = LNAsim_sensitivities_sge(paramsSGE(1:4));
[x2] = LNAsim_sensitivities_sge(paramsNF1000(1:4));
[x3] = LNAsim_sensitivities_sge(paramsNF200(1:4));
[x4] = LNAsim_sensitivities_sge(paramsNF50(1:4));

% rna level
observables(:,1) = [x(1),x2(1),x3(1),x4(1)];
% protein level
observables(:,2) = [x(2),x2(2),x3(2),x4(2)];
% RNA variance
observables(:,3) = [x(3),x2(3),x3(3),x4(3)];
% protein variance
observables(:,4) = [x(4),x2(4),x3(4),x4(4)];
% covariance
observables(:,5) = [x(5),x2(5),x3(5),x4(5)];
% RNA CV
observables(:,6) = [sqrt(x(3))/x(1),sqrt(x2(3))/x2(1),sqrt(x3(3))/x3(1),sqrt(x4(3))/x4(1)];
% Protein CV
observables(:,7) = [sqrt(x(4))/x(2),sqrt(x2(4))/x2(2),sqrt(x3(4))/x3(2),sqrt(x4(4))/x4(2)];

all_observables{2} = observables;

%% Second: parameters effect on NF with hill 2

[x] = LNAsim_sensitivities_sge_nfeedback_hill2(paramsSGE);
[x2] = LNAsim_sensitivities_sge_nfeedback_hill2(paramsNF1000);
[x3] = LNAsim_sensitivities_sge_nfeedback_hill2(paramsNF200);
[x4] = LNAsim_sensitivities_sge_nfeedback_hill2(paramsNF50);

observables = zeros(4,7);

% rna level
observables(:,1) = [x(1),x2(1),x3(1),x4(1)];
% protein level
observables(:,2) = [x(2),x2(2),x3(2),x4(2)];
% RNA variance
observables(:,3) = [x(3),x2(3),x3(3),x4(3)];
% protein variance
observables(:,4) = [x(4),x2(4),x3(4),x4(4)];
% covariance
observables(:,5) = [x(5),x2(5),x3(5),x4(5)];
% RNA CV
observables(:,6) = [sqrt(x(3))/x(1),sqrt(x2(3))/x2(1),sqrt(x3(3))/x3(1),sqrt(x4(3))/x4(1)];
% Protein CV
observables(:,7) = [sqrt(x(4))/x(2),sqrt(x2(4))/x2(2),sqrt(x3(4))/x3(2),sqrt(x4(4))/x4(2)];

all_observables{3} = observables;


%% 
% Second: optimisation on covariance

observables = zeros(4,7);
params2 = [10,0.004,0.001,0.0025,200];
params4 = [10,0.0035,0.001,0.0028,200];
params6 = [10,0.0034,0.001,0.0033,200];
params8 = [10,0.0033,0.001,0.0031,200];

[x2] = LNAsim_sensitivities_sge_nfeedback_hill2(params2);
[x4] = LNAsim_sensitivities_sge_nfeedback_hill4(params4);
[x6] = LNAsim_sensitivities_sge_nfeedback_hill6(params6);
[x8] = LNAsim_sensitivities_sge_nfeedback_hill8(params8);

% rna level
observables(:,1) = [x2(1),x4(1),x6(1),x8(1)];
% protein level
observables(:,2) = [x2(2),x4(2),x6(2),x8(2)];
% RNA variance
observables(:,3) = [x2(3),x4(3),x6(3),x8(3)];
% protein variance
observables(:,4) = [x2(4),x4(4),x6(4),x8(4)];
% covariance
observables(:,5) = [x2(5),x4(5),x6(5),x8(5)];
% RNA CV
observables(:,6) = [sqrt(x2(3))/x2(1),sqrt(x4(3))/x4(1),sqrt(x6(3))/x6(1),sqrt(x8(3))/x8(1)];
% Protein CV
observables(:,7) = [sqrt(x2(4))/x2(2),sqrt(x4(4))/x4(2),sqrt(x6(4))/x6(2),sqrt(x8(4))/x8(2)];

all_observables{3} = observables;
