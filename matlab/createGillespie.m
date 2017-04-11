%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                Generalised function to create Gillespie             %
%                   simulation file from input file                   %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Takes as input the circuit name, the indices of any dimerisation reaction
%%% there may be, as well as the index of the species involved in the 
%%% dimers - e.g. createGillespie('pfeedback',[5],[3]) for dimerisation
%%% of species 3 by reaction 5.

function createGillespie(circuitName,dimerReactions,dimerSpecies)

    str1 = 'input_';
    filename = strcat(str1,circuitName);

    %Call file 'input_circuitName' which has all rate constants and species
    %organised into vector 'const' and matrix 'species'

    eval(filename);

    nspecies = length(xinit);
    nreactions = length(macroF);
    nparam = length(param);
    
    %% Create a vector with symbolic variable names for each species
  
    str1a = 'x';
    symX = sym(zeros(1,nspecies));
  
    for i = 1:nspecies   
    
        str1b = int2str(i);
        symX(i) = sym(strcat(str1a,str1b,str1a));
    
    end
    
    %% Convert rates vector macroF into numeric variables
  
    % Substitute 'x1x' type variables by 'x(1)' 
  
    macroFnew = macroF;
    
    for k = 1:length(symX)
      
      oldvar = strcat('x',num2str(k),'x');
      newvar = strcat('species(',num2str(k),')');
      
      for p = 1:length(macroF)
          
          oldrate = char(macroFnew(p));
          newrate = strrep(oldrate,oldvar,newvar);
          macroFnew(p) = sym(newrate);
          
      end
      
    end
    
    %% Create the .m file in a convenient layout
  
    cellstring = {};
    
    cellstring(1) = {strcat('function [timepoints,species_timeseries] = gillespie_',circuitName,'(simLength)')};
    cellstring(2) = {''};
    cellstring(3) = {strcat('initialConditions = [',num2str(xinit),'];')};
    cellstring(4) = {strcat('nspecies = length(initialConditions);')};
    cellstring(5) = {strcat('nreactions = ',num2str(nreactions),';')};
    cellstring(6) = {'species = initialConditions;'};
    cellstring(7) = {'timepoints = zeros(1000000,1);'};
    cellstring(8) = {'species_timeseries = zeros(length(timepoints),nspecies);'};
    cellstring(9) = {''};
    cellstring(10) = {'S = zeros(nspecies,nreactions);'};
    
    cellLength = length(cellstring);
    
    for i = 1:nspecies
        
        cellstring(cellLength + i) = {strcat('S(',num2str(i),',:) = [',num2str(S(i,:)),'];')};
            
    end
    
    cellLength = length(cellstring);
    
    cellstring(cellLength + 1) = {''};
    cellstring(cellLength + 2) = {'%% Input parameters'};
    cellstring(cellLength + 3) = {''};
    
    cellLength = length(cellstring);
    
    for i = 1:nparam
        
        cellstring(cellLength + i) = {strcat(char(param(i)),' = ',num2str(eval(char(param(i)))),';')};
    
    end
    
    cellLength = length(cellstring);
    
    cellstring(cellLength + 1) = {''};
    cellstring(cellLength + 2) = {'ttot = 0;'};
    cellstring(cellLength + 3) = {''};
    cellstring(cellLength + 4) = {'% Initialise reaction rates'};
    
    cellLength = length(cellstring);
    
    for i = 1:nreactions
    
        cellstring(cellLength + i) = {strcat('rates(',num2str(i),') = ',char(macroFnew(i)),';')};
        
    end
    
    cellLength = length(cellstring);    
    
    count = 0;
    
    for j = 1:length(dimerReactions)
    
    count = count + 3;
    cellstring(cellLength + count) = {strcat('if(species(',num2str(dimerSpecies(j)),') < 2)')};
    cellstring(cellLength + count + 1) = {strcat('rates(',num2str(dimerReactions(j)),') = 0;')};
    cellstring(cellLength + count + 2) = {'end'};
    
    end

    cellLength = length(cellstring);

    cellstring(cellLength + 1) = {''};
    cellstring(cellLength + 2) = {'counter = 0;'};
    cellstring(cellLength + 3) = {'breaktest = 0;'};
    cellstring(cellLength + 4) = {''};
    cellstring(cellLength + 5) = {'while(ttot<simLength)'};
    cellstring(cellLength + 6) = {''};
    cellstring(cellLength + 7) = {'counter = counter + 1;'};
    cellstring(cellLength + 8) = {''};
    cellstring(cellLength + 9) = {'mrands = rand(1,nreactions);'};
    cellstring(cellLength + 10) = {''};
    cellstring(cellLength + 11) = {'theta = (1./rates).*log(1./mrands);'};
    cellstring(cellLength + 12) = {'negtest = find(theta==-Inf);'};
    cellstring(cellLength + 13) = {'theta(negtest) = Inf;'};
    cellstring(cellLength + 14) = {''};
    cellstring(cellLength + 15) = {'tfirst = min(theta);'};
    cellstring(cellLength + 16) = {'ttot = ttot + tfirst;'};
    cellstring(cellLength + 17) = {'if(ttot==-Inf)'};
    cellstring(cellLength + 18) = {'breaktest = 1;'};
    cellstring(cellLength + 19) = {'break'};
    cellstring(cellLength + 20) = {'end'};
    cellstring(cellLength + 21) = {'timepoints(counter) = ttot;'};
    cellstring(cellLength + 22) = {'index = find(theta==tfirst);'};
    cellstring(cellLength + 23) = {''};
    cellstring(cellLength + 24) = {'species = species + transpose(S(:,index));'};
    cellstring(cellLength + 25) = {''};
    cellstring(cellLength + 26) = {'% Update rates'};
    
    cellLength = length(cellstring);
    
    for i = 1:nreactions
    
        cellstring(cellLength + i) = {strcat('rates(',num2str(i),') = ',char(macroFnew(i)),';')};
        
    end
    
    cellLength = length(cellstring);
    
    count = 0;
    
    for j = 1:length(dimerReactions)
    
    count = count + 3;
    cellstring(cellLength + count) = {strcat('if(species(',num2str(dimerSpecies(j)),') < 2)')};
    cellstring(cellLength + count + 1) = {strcat('rates(',num2str(dimerReactions(j)),') = 0;')};
    cellstring(cellLength + count + 2) = {'end'};
    
    end

    cellLength = length(cellstring);
    
    cellstring(cellLength + 1) = {''};
    cellstring(cellLength + 2) = {'% Append species to species_timeseries'};
    cellstring(cellLength + 3) = {''};
    cellstring(cellLength + 4) = {'species_timeseries(counter,:) = species;'};
    cellstring(cellLength + 5) = {''};
    cellstring(cellLength + 6) = {'end'};
    cellstring(cellLength + 7) = {''};
    cellstring(cellLength + 8) = {'% Truncate all zeros at the end of the timeseries'};
    cellstring(cellLength + 9) = {'species_timeseries = species_timeseries(1:counter,:);'};
    cellstring(cellLength + 10) = {'timepoints = timepoints(1:counter);'};
    
    cellLength = length(cellstring);
    
    str1 = str2mat('touch gillespie_');
    str2 = str2mat(circuitName);
    str3 = str2mat('.m');

    sysCommand = [str1,str2,str3];
    
    system(sysCommand);
  
    filename = strcat('gillespie_',circuitName,'.m');
  
    fid = fopen(filename,'w');
    
    for i = 1:cellLength 
        
        temp = cell2mat(cellstring(i));
        fprintf(fid,[temp '\n']);
        
    end
        
    fclose(fid);
    
end


    
    
