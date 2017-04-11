function initConditionsCell = generateICs_p1000(x0,nconditions)

k1Vals = logspace(0,1,nconditions);
k2Vals = logspace(-1,1,nconditions);
k4Vals = [0.01 0.001];

paramVals = zeros(nconditions,length(x0));
counter = 0;

for i = 1:length(k1Vals)
    
    for j = 1:length(k2Vals)
        
        for k = 1:length(k4Vals)
            
            counter = counter + 1;
            newparams = [k1Vals(i),k2Vals(j),1000*k2Vals(j)*k4Vals(k)/k1Vals(i),k4Vals(k)];
            paramVals(counter,:) = newparams;
            
        end
        
    end
    
end

dim = size(paramVals);
numOfInitConditions = dim(1);

initConditionsCell = {};

for j = 1:numOfInitConditions
    
    initConditionsCell{j} = {paramVals(j,:)};
    
end

end