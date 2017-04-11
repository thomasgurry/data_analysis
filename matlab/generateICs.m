%% Generate many different initial conditions from 10% to 1000% of original
%% value.  Returns a cell array.  Note that the input parameter 'nconditions' 
%% specifies the number of different conditions PER PARAMATER (e.g. 3 
%% parameters and nconditions = 2 results in 2^3 combinations).

function initConditionsCell = generateICs(x0,nconditions)

numOfInitConditions = nconditions^length(x0);
initConditions = zeros(numOfInitConditions,length(x0));
paramRange = logspace(-1,1,nconditions);
paramVals = zeros(nconditions,length(x0)); % stores all possible values of params

for i = 1:length(x0)
   
    currentBaseline = x0(i);
    paramVals(:,i) = paramRange*currentBaseline;
    
end

% Generate matrix of initial conditions
    
commandStr = 'initConditions = setprod(';

for i = 1:(length(x0)-1)
    
    commandStr = strcat(commandStr,'paramVals(:,',num2str(i),'),');
    
end

commandStr = strcat(commandStr,'paramVals(:,',num2str(length(x0)),'));');

eval(commandStr)

initConditionsCell = {};

for j = 1:numOfInitConditions
    
    initConditionsCell{j} = {initConditions(j,:)};
    
end

end