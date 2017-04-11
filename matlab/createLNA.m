%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                Generalised function to create ODE file to           %
%                   simulate a system using the Linear                %
%                          Noise Approximation                        %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Takes as input the name of the circuit as a string, and the name of the
%%% M file with parameters, rates and stoichiometry matrix, also as a string
%%% (e.g. 'input_sge').

function createLNA(circuitName)

  str1 = 'input_';
  inputFile = strcat(str1,circuitName);
  eval(inputFile)  

  nspecies = length(xinit);
  nrates = length(macroF);
  nElements = nspecies^2;  %number of elements in A, EEt and C
  
  %% Create a vector with symbolic variable names for each species
  
  str1a = 'x';
  symX = sym(zeros(1,nspecies));
  
  for i = 1:nspecies   
    
    str1b = int2str(i);
    symX(i) = sym(strcat(str1a,str1b,str1a));
    
  end
  
  %% Compute the Jacobian matrix J of the rates wrt the species 
  %% as a symbolic matrix
  
  J = sym(zeros(nrates,nspecies));
  
  for j = 1:nrates
    
    for k = 1:nspecies
    
      J(j,k) = diff(macroF(j),symX(k));
      
    end
    
  end
      
      
  %% Compute matrices A, E and E*transpose(E) in symbolic forms
  
  symS = sym(S);
  symA = symS*J;  
  symE = symS*sqrt(diag(macroF));
  symEEt = symE*transpose(symE);
  
  %% Compute the covariance matrix C in symbolic form - c1c, c2c etc.
  
  str3 = 'c';
  symC = sym(zeros(nspecies));
  
  for i = 1:nElements   
    
    str4 = int2str(i);
    symC(i) = sym(strcat(str3,str4,str3));
    
  end
  
  % Make symC symmetric
  tempMat = triu(symC);
  tempMatdiag = diag(diag(symC));
  tempMat = tempMat + transpose(tempMat) - tempMatdiag;
  symC = tempMat;
  
  % Set up dC/dt = AC + CA^T + (EE^T) 
  mat1 = symA*symC + symC*transpose(symA) + symEEt;
  
  covMat_dot = mat1; 
  
  for l = 1:length(symX)
      
      oldvar = strcat('x',num2str(l),'x');
      newvar = strcat('x(',num2str(l),')');
      
      for q = 1:nElements
          
          oldcov = char(covMat_dot(q));
          newcov = strrep(oldcov,oldvar,newvar);
          covMat_dot(q) = sym(newcov);
          
      end
      
  end
  
  % Create initial covariance matrix by setting the standard deviations of
  % the initial species numbers to be 10% of their value.  Covariances set
  % to 0.  Corresponds to an experiment in which there is 10% standard
  % deviation in the initial conditions (e.g. pipetting
  % error).  Furthermore, allows for 0 variance when species number is 0.
  
  Czero = zeros(size(covMat_dot));
  
  for i = 1:nspecies
    
    Czero(i,i) = (xinit(i)/10)^2;
    
  end
  
  %% Substitute c1, c2, c3 with x(nspecies+1), etc. in covMat_dot
  
  % First the variances
  
  for i = 1:nspecies
    
    oldvar = char(symC(i,i));
    newvar = strcat('x(',num2str(nspecies+i),')');
    
    for j = 1:nElements
  
      oldcov = char(covMat_dot(j));
      newcov = strrep(oldcov,oldvar,newvar);
      covMat_dot(j) = sym(newcov);
      
    end
    
  end
  
  % Then the covariances
  
  counter = 0;
  
  for i = 1:nspecies
    
    for j = (i+1):nspecies
      
      counter = counter + 1;
      oldvar = char(symC(i,j));
      newvar = strcat('x(',num2str(nspecies+nspecies+counter),')');
      
      for k = 1:nElements
	
	oldcov = char(covMat_dot(k));
	newcov = strrep(oldcov,oldvar,newvar);
	covMat_dot(k) = sym(newcov);
	
      end
      
    end
    
  end
  	 
  
  %% Append initial variances and covariances to vector 'xinit' for input
  %% into ODE solver

  counter = length(xinit);
  
  % First the variances
  
  for i = 1:nspecies
    
    counter = counter + 1;
    xinit(counter) = Czero(i,i);
    
  end
  
  % Then the covariances
  
  for i = 1:nspecies
    
    for j = (i+1):nspecies
      
      counter = counter + 1;
      xinit(counter) = Czero(i,j);
      
    end
    
  end			   
  
  %% Convert rates vector macroF into numeric variables with parameters
  %% substituted
  
   % Substitute all parameters in 'macroF' by their numeric values
  
  macroFnew = macroF;
   
  for i = 1:length(param)
    
    oldvar = char(param(i));
    newvar = strcat(num2str(eval(oldvar)));
    
    for q = 1:length(macroFnew)
      
      oldrate = char(macroFnew(q));
      newrate = strrep(oldrate,oldvar,newvar);
      macroFnew(q) = sym(newrate);
      
    end
    
  end
  
  % Substitute 'x1x' type variables by 'x(1)' 
  
  for k = 1:length(symX)
      
      oldvar = strcat('x',num2str(k),'x');
      newvar = strcat('x(',num2str(k),')');
      
      for p = 1:length(macroF)
          
          oldrate = char(macroFnew(p));
          newrate = strrep(oldrate,oldvar,newvar);
          macroFnew(p) = sym(newrate);
          
      end
      
  end
  
  %% Substitute all parameters in covMat_dot by their numeric values
  
  covMat_dot_new = covMat_dot;
   
  for i = 1:length(param)
    
    oldvar = char(param(i));
    newvar = strcat(num2str(eval(oldvar)));
    
    for q = 1:nElements
      
      oldcov = char(covMat_dot_new(q));
      newcov = strrep(oldcov,oldvar,newvar);
      covMat_dot_new(q) = sym(newcov);
      
    end
    
  end
  
  covMat_dot = covMat_dot_new;
  
  %% Create a file 'MRE_circuitName.m' to be used by the ODE solver.
  %% Contains S, macroF and covMat_dot. 
  
  str1 = str2mat('touch MRE_');
  str2 = str2mat(circuitName);
  str3 = str2mat('.m');

  sysCommand = [str1,str2,str3];
    
  system(sysCommand);
  
  filename = strcat('MRE_',circuitName,'.m');
  
  fid = fopen(filename,'w');
  
  cellstring = {};
  
  cellstring{1} = strcat('function dxdt = MRE_',circuitName,'(t,x)');
  cellstring{2} = '';
  cellstring{3} = strcat('nspecies = ',num2str(nspecies),';');
  cellstring{4} = strcat('S = zeros(nspecies,',num2str(nrates),');');
  
  nElementsS = nspecies*nrates;
  
  % Write S to file
  
  for k = 1:nElementsS
      
      cellstring{4+k} = strcat('S(',num2str(k),') = ',num2str(S(k)),';');
      
  end
  
  cellLength = length(cellstring);
  
  cellstring{cellLength+1} = strcat('rates = zeros(1,',num2str(nrates),');');
  
  cellLength = length(cellstring);
  
  % Write macroFnew to file
  
  for j = 1:length(macroF)
      
      cellstring{cellLength+j} = strcat('rates(',num2str(j),') = ',char(macroFnew(j)),';');
      
  end
  
  
  % Write covMat_dot to file, as a vector with the variances first,
  % followed by the covariances (only upper triangle), row by row.
  
  nVarAndCovar = (nspecies^2 - nspecies)/2 + nspecies;
  cellstring{length(cellstring) + 1} = strcat('covMat_dot = zeros(1,',num2str(nVarAndCovar),');');
  cellLength = length(cellstring);
  
  for i = 1:nspecies
    
    cellstring{cellLength + i} = strcat('covMat_dot(',num2str(i),') = ',char(covMat_dot(i,i)),';');   
  
  end
  
  cellLength = length(cellstring);
  
  counter = 0;
  
  for i = 1:nspecies
  
    for j = i+1:nspecies
      
      counter = counter + 1;
    
      cellstring{cellLength + counter} = strcat('covMat_dot(',num2str(nspecies + counter),') = ',char(covMat_dot(i,j)),';');
      
    end
    
  end
  
  cellLength = length(cellstring);
  
  cellstring{cellLength + 1} = '';
  
  % 'nterms' is the number of terms in the ODE, the first length(x) are the
  % species for the MRE, and the next (length(x)^2 - length(x))/2 + length(x) are the
  % variances and covariances
  
  cellstring{cellLength + 2} = 'nterms = length(x);';
  cellstring{cellLength + 3} = 'dxdt = zeros(nterms,1);';
  cellstring{cellLength + 4} = '';
  cellstring{cellLength + 5} = 'for i = 1:nspecies';
  cellstring{cellLength + 6} = '';
  cellstring{cellLength + 7} = 'dxdt(i) = S(i,:)*transpose(rates);';
  cellstring{cellLength + 8} = '';
  cellstring{cellLength + 9} = 'end';
  cellstring{cellLength + 10} = '';
  cellstring{cellLength + 11} = 'counter = 0;';
  cellstring{cellLength + 12} = '';
  cellstring{cellLength + 13} = 'for j = (nspecies+1):nterms';
  cellstring{cellLength + 14} = '';
  cellstring{cellLength + 15} = 'counter = counter + 1;';
  cellstring{cellLength + 16} = 'dxdt(j) = covMat_dot(counter);';
  cellstring{cellLength + 17} = '';
  cellstring{cellLength + 18} = 'end';
  
  for i = 1:length(cellstring)
      
      temp = cell2mat(cellstring(1,i));
      fprintf(fid,[temp '\n']);
      
  end
  
  fclose(fid);
  
  
%% Create a file LNAsim_.m which has the initial conditions of each
%% variable (all the species, variances and covariances) and the commands
%% to simulate the system

    str1 = str2mat('touch LNAsim_');
    str2 = str2mat(circuitName);
    str3 = str2mat('.m');

    sysCommand = [str1,str2,str3];

    system(sysCommand);

    filename = strcat('LNAsim_',circuitName,'.m');

    fid = fopen(filename,'w');

    cellstring = {};

    cellstring{1} = strcat('function [timepoints,Xvals] = LNAsim_',circuitName,'(simLength)');
    cellstring{2} = '';

    for i = 1:length(xinit)

        cellstring{2+i} = strcat('xinit(',num2str(i),') = ',num2str(xinit(i)),';');

    end
    
    cellLength = length(cellstring);
    
    functionName = strcat('MRE_',circuitName);
    
    cellstring{cellLength + 1} = 'n = simLength*100;';
    cellstring{cellLength + 2} = 'tspan = linspace(0,simLength,n);';
    cellstring{cellLength + 3} = strcat('[timepoints,Xvals] = ode15s(@',functionName,',tspan,xinit);');
    
    for i = 1:length(cellstring)
        
        temp = cell2mat(cellstring(1,i));
        fprintf(fid,[temp '\n']);
        
    end
    
    fclose(fid);
    
end
