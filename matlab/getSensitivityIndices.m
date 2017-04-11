%% Script to calculate which indices represent the first and last indices
%% of the parameter sensitivities of a given species.

%% Takes as input the number of species (NOT including variances etc)
%% 'nspecies', the number of parameters 'nparams' and the species index
%% 'sindex'

function [inds] = getSensitivityIndices(nspecies,nparams,sindex)

n = nspecies;
p = nparams;
q = sindex;

inds = [0,0];
inds(1) = n + n*(n+1)/2 + (q-1)*p + 1;
inds(2) = n + n*(n+1)/2 + q*p;

end