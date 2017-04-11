%% Function to append the equilibrated sub-system to the
%% whole system's initial conditions

function [yinit] = LNAappend_init_sensitivities(x,nspecies_subsystem,nparams_subsystem,nspecies_fullsystem,nparams_fullsystem)

% Add all equilibrated sub-system to full system initial conditions

[species,covmat,sensitivityMat] = LNAvect2mat(x,nspecies_subsystem,nparams_subsystem);
species_full = zeros(1,nspecies_fullsystem);
covmat_full = zeros(nspecies_fullsystem,nspecies_fullsystem);
nVarsAndCovars_fullsystem = nspecies_fullsystem*(nspecies_fullsystem+1)/2;
sensitivityMat_full = zeros(nparams_fullsystem,nspecies_fullsystem + nVarsAndCovars_fullsystem);

species_full(1:nspecies_subsystem) = species;
covmat_full(1:nspecies_subsystem,1:nspecies_subsystem) = covmat;
sensitivityMat_full(1:nparams_subsystem,1:size(sensitivityMat,2)) = sensitivityMat;

% Convert into a single vector yinit

yinit = [];
n = nspecies_fullsystem; % number of species
N = nspecies_fullsystem + nVarsAndCovars_fullsystem; % number of species + their vars & covars
p = nparams_fullsystem;
yinit(1:n) = species_full;
counter = n;

for i = 1:n
    
    counter = counter + 1;
    yinit(counter) = covmat_full(i,i);  % variances
    
end

for i = 1:n
    
    for j = (i+1):n
        
        counter = counter + 1;
        yinit(counter) = covmat_full(i,j);  % covariances
        
    end
    
end

% Sensitivities sorted species by species, i.e. [dx1dp,dx2dp,dx3dp,etc.]

for i = 1:p
    
    for j = 1:N
        
        counter = counter + 1;
        yinit(counter) = sensitivityMat_full(i,j);
        
    end
    
end

end


%% Function that takes as input the vector of species, variances,
%% covariances and sensitivities that is the output of an LNA simulation
%% and converts it vector of species, full covariance matrix and full 
%% sensitivity matrix.

function [species,covmat,sensitivityMat] = LNAvect2mat(x,nspecies,nparams)

species = x(1:nspecies);
covmat = zeros(nspecies);
nVarsAndCovars = nspecies*(nspecies+1)/2;
sensitivityMat = zeros(nparams,nspecies + nVarsAndCovars); 

% Get indices of variances and covariances

varInds_vect = [nspecies+1 : 2*nspecies];
covarInds_vect = [2*nspecies + 1 : nspecies + nspecies*(nspecies+1)/2];

variances_mat = diag(x(varInds_vect));
covariances_mat = zeros(nspecies);
counter = 0;

for i = 1:nspecies
    
    for j = (i+1):nspecies
        
        counter = counter + 1;
        covariances_mat(i,j) = x(covarInds_vect(counter));
        
    end
    
end

covmat = variances_mat + covariances_mat + covariances_mat.';

% Fill sensitivity matrix

sensitivityInds_vect = [nspecies + nspecies*(nspecies+1)/2 + 1 : length(x)];
sensitivityMat(:) = x(sensitivityInds_vect);

end