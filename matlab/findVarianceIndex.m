%% Function to give you the index of the variance of a given species
%% from a LNA Kronecker simulation 'sim' using the 'SimulateLna' function

function index = findVarianceIndex(speciesIndex, nspecies)

index = nspecies;
n = speciesIndex;

% recursively find index of the nth species variance

for i = 1:n
    
    index = index + i;
    
end

end