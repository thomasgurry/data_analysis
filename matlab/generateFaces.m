%% Script to compute the parameter values for +/- 10% of each
%% parameter value individually and output the set of new parameters.

function paramCombinations = generateFaces(x0)

numofParamCombinations = 2*length(x0);
paramCombinations = zeros(numofParamCombinations,length(x0));

% Generate +10% and -10% values for each parameter

faceValues = zeros(2,length(x0));
faceValues(1,:) = 1.1*x0;
faceValues(2,:) = 0.9*x0;

% Generate matrix of face values

counter = 0;

for i = 1:length(x0)
    
    counter = counter + 1;
    paramCombinations(counter,:) = x0;
    paramCombinations(counter,i) = faceValues(1,i);
    
    counter = counter + 1;
    paramCombinations(counter,:) = x0;
    paramCombinations(counter,i) = faceValues(2,i);
    
end

end