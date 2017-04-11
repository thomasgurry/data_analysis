function [] = ensemble_RG_distribution(rgVals, weights)

[N,binIndices] = histc(rgVals,[1:70]);

binWeights = zeros(1,70);

for i = 1:70
    
    binInds = find(binIndices==i);
    binWeights(i) = sum(weights(binInds));
    
end

% Ensure rgVals is row vector and weights a column vector
if(size(rgVals,1)>1)
    rgVals = rgVals.';
end

if(size(weights,2)>1)
    weights = weights.';
end

ensemble_avg_RG = rgVals*weights;

figure
subplot(2,1,1)
plot([1:70],binWeights,'-ok')
xlabel('Rg (Angstroms)')
ylabel('Probability density')
title('Ensemble distribution of radius of gyration')
hold on
line([ensemble_avg_RG,ensemble_avg_RG],[0,max(binWeights)])
hold off

subplot(2,1,2)
plot([1:70],smooth(binWeights),'-ok')
xlabel('Rg (Angstroms)')
ylabel('Probability density')
title('Ensemble distribution of radius of gyration (smoothed)')
hold on
line([ensemble_avg_RG,ensemble_avg_RG],[0,max(binWeights)])
hold off

end