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
tic
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
  
  % Store covMat_dot in its current form as 'originalMat' for use later in calculating the jacobian of the LNA
  % system.
  
  originalMat = covMat_dot;
  
  for l = 1:length(symX)
      
      oldvar = strcat('x',num2str(l),'x');
      newvar = strcat('x(',num2str(l),')');
      
      for q = 1:nElements
          
          oldcov = char(covMat_dot(q));
          newcov = strrep(oldcov,oldvar,newvar);
          covMat_dot(q) = sym(newcov);
          
      end
      
  end
			 
  toc
  
  % Create initial covariance matrix by setting the standard deviations of
  % the initial species numbers to be 10% of their value.  Covariances set
  % to 0.  Corresponds to an experiment in which there is 10% standard
  % deviation in the administration of initial conditions (e.g. pipetting
  % error).  Furthermore, allows for 0 variance when species number is 0.
  
  Czero = zeros(size(covMat_dot));
  
  for i = 1:nspecies
    
    Czero(i,i) = (xinit(i)/10)^2;
    
  end
  
  %% Substitute c1, c2, c3 with x(nspecies+1), etc. in covMat_dot and
  %% originalMat
  
  % First the variances
  
  for i = 1:nspecies
    
    oldvar = char(symC(i,i));
    newvar1 = strcat('x(',num2str(nspecies+i),')');
    newvar2 = strcat('x',num2str(nspecies+i),'x');
    
    for j = 1:nElements
  
      oldcov1 = char(covMat_dot(j));
      newcov1 = strrep(oldcov1,oldvar,newvar1);
      covMat_dot(j) = sym(newcov1);
      
      oldcov2 = char(originalMat(j));
      newcov2 = strrep(oldcov2,oldvar,newvar2);
      originalMat(j) = sym(newcov2);
      
    end
    
  end
  
  % Then the covariances
  
  counter = 0;
  
  for i = 1:nspecies
    
    for j = (i+1):nspecies
      
      counter = counter + 1;
      
      oldvar = char(symC(i,j));
      newvar1 = strcat('x(',num2str(nspecies+nspecies+counter),')');
      newvar2 = strcat('x',num2str(nspecies+nspecies+counter),'x');
      
      for k = 1:nElements
	
        oldcov1 = char(covMat_dot(k));
        newcov1 = strrep(oldcov1,oldvar,newvar1);
        covMat_dot(k) = sym(newcov1);
    
        oldcov2 = char(originalMat(k));
        newcov2 = strrep(oldcov2,oldvar,newvar2);
        originalMat(k) = sym(newcov2);
	
      end
      
    end
    
  end
  
  % Update the vector 'symX' containing the symbolic name of all species.
  
  nspeciesTotal = nspecies+nspecies+counter;
  
  symX = sym(zeros(1,nspeciesTotal));
  
  for i = 1:nspeciesTotal
      
      symX(i) = strcat('x',num2str(i),'x');
      
  end
  
  toc
  
  %% Store all LNA equations in one array 'LNAequations' (Master Rate Equations + 
  %% Covariance Matrix Equations), and compute the Jacobian of the LNA 
  %% equations with prespect to parameters and with respect to species 
  
  LNAequations = sym(zeros(1));
  
  % Master rate equations
  
  for i = 1:nspecies
      
      LNAequations(i) = symS(i,:)*transpose(macroF);
      
  end
  
  % Variance equations
  
  counter = nspecies;
  
  for i = 1:nspecies
      
    counter = counter + 1;
    LNAequations(counter) = originalMat(i,i);
          
  end
  
  % Covariance equations
  
  for i = 1:nspecies 
      
      for j = (i+1):nspecies
          
          counter = counter + 1;
          LNAequations(counter) = originalMat(i,j);
          
      end
      
  end
 
  % Jacobian with respect to param
  
  param_jacobian = sym(zeros(length(LNAequations),length(param)));
  
  for i = 1:length(LNAequations)
      
      for j = 1:length(param)
  
        param_jacobian(i,j) = diff(LNAequations(i),char(param(j)));
        
      end
      
  end
  
  % Jacobian with respect to species
  
  species_jacobian = sym(zeros(length(LNAequations),nspecies));
  
  for i = 1:length(LNAequations)
      
      for j = 1:length(symX)
          
          species_jacobian(i,j) = diff(LNAequations(i),char(symX(j)));
          
      end
      
  end
  
  %% Obtain the system of ODE's for the sensitivites.  If F = dx/dt, then
  %% we have d(dx/dp)/dt = dF/dp + (dF/dx)*(dx/dp), or dY/dt = A + BY,
  %% where A is the Jacobian of the LNA equations with respect to
  %% parameters and B the Jacobian of the LNA equations with respect to
  %% species.
  
  toc
  
  dimSensitivities = size(param_jacobian);
  nElementsYdot = length(LNAequations)*length(param);
  Y_dot = sym(zeros(dimSensitivities));
  Y_mat = sym(zeros(dimSensitivities));
  counter = 0;
  str1 = 'y';
  
  for i = 1:dimSensitivities(1)
      
      for j = 1:dimSensitivities(2)
          
          counter = counter + 1;
          Y_mat(i,j) = strcat(str1,num2str(counter),str1);
          
      end
      
  end
  
  Y_dot = param_jacobian + species_jacobian*Y_mat;
  
  toc
  
  % Replace all terms 'x1x' by 'x(1)' in Y_dot
  for i = 1:length(symX)
      
      oldvar = char(symX(i));
      newvar = strcat('x(',num2str(i),')');
      
      for j = 1:nElementsYdot
          
          oldterm = char(Y_dot(j));
          newterm = strrep(oldterm,oldvar,newvar);
          Y_dot(j) = sym(newterm);
          
      end
      
  end
  
  toc
  
  % Replace all parameters 'k1' by 'y(1)'
  for i = 1:length(param)
      
      oldvar = char(param(i));
      newvar = strcat('y(',num2str(i),')');
      
      for j = 1:nElementsYdot
          
          oldterm = char(Y_dot(j));
          newterm = strrep(oldterm,oldvar,newvar);
          Y_dot(j) = sym(newterm);
          
      end
      
  end
  
  
  toc
  
  % Replace all 'y1y' (parameter sensitivities) with 'x(N+1)'
  N = length(symX);
  counter = 0;
  
   for i = 1:nElementsYdot
      
      counter = counter + 1;
      oldvar = strcat('y',num2str(counter),'y');
      newvar = strcat('x(',num2str(N+counter),')');
      
      for j = 1:nElementsYdot
          
          oldterm = char(Y_dot(j));
          newterm = strrep(oldterm,oldvar,newvar);
          Y_dot(j) = sym(newterm);
          
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
  
  %% Substitute parameters by 'k1' by 'y(1)' in rates 'macroFnew' and 
  %% by their numeric value in 'macroFparam' (numeric values specified in
  %% the input file).
  
  toc
  
  macroFnew = macroF;
  macroFparam = macroF;
   
  for i = 1:length(param)
    
    oldvar = char(param(i));
    newvar1 = strcat('y(',num2str(i),')');
    newvar2 = strcat(num2str(eval(oldvar)));
    
    for q = 1:length(macroF)
      
      oldrate1 = char(macroFnew(q));
      newrate1 = strrep(oldrate1,oldvar,newvar1);
      macroFnew(q) = sym(newrate1);
      
      oldrate2 = char(macroFparam(q));
      newrate2 = strrep(oldrate2,oldvar,newvar2);
      macroFparam(q) = sym(newrate2);
      
    end
    
  end
  
  toc
  
  % Substitute 'x1x' type variables by 'x(1)' 
  
  for k = 1:length(symX)
      
      oldvar = strcat('x',num2str(k),'x');
      newvar = strcat('x(',num2str(k),')');
      
      for p = 1:length(macroF)
          
          oldrate1 = char(macroFnew(p));
          oldrate2 = char(macroFparam(p));
          newrate1 = strrep(oldrate1,oldvar,newvar);
          newrate2 = strrep(oldrate2,oldvar,newvar);
          macroFnew(p) = sym(newrate1);
          macroFparam(p) = sym(newrate2);
          
      end
      
  end
  
  %% Substitute all parameters in covMat_dot to 'y(1)' type and create a
  %% covMat_dot_param containing the numeric values of the parameters.
  
  covMat_dot_new = covMat_dot;
  covMat_dot_param = covMat_dot;
   
  for i = 1:length(param)
    
    oldvar = char(param(i));
    newvar1 = strcat('y(',num2str(i),')');
    newvar2 = strcat(num2str(eval(oldvar)));
    
    for q = 1:nElements
      
      oldcov1 = char(covMat_dot_new(q));
      newcov1 = strrep(oldcov1,oldvar,newvar1);
      covMat_dot_new(q) = sym(newcov1);
      
      oldcov2 = char(covMat_dot_param(q));
      newcov2 = strrep(oldcov2,oldvar,newvar2);
      covMat_dot_param(q) = sym(newcov2);
      
    end
    
  end
  
  covMat_dot = covMat_dot_new;
  
  toc
  
  %% Append values of 0 to initial parameter sensitivities in 'xinit'
  
  xinit_no_sensitivities = xinit;
  xinit_with_sensitivities = xinit;
  lengthXinit = length(xinit);
  
  for i = 1:nElementsYdot
      
      xinit_with_sensitivities(lengthXinit + i) = 0;
      
  end
  
  
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
      
      cellstring{cellLength+j} = strcat('rates(',num2str(j),') = ',char(macroFparam(j)),';');
      
  end
  
  % Write covMat_dot to file, as a vector with the variances first,
  % followed by the covariances (only upper triangle), row by row.
  
  nVarAndCovar = (nspecies^2 - nspecies)/2 + nspecies;
  cellstring{length(cellstring) + 1} = strcat('covMat_dot = zeros(1,',num2str(nVarAndCovar),');');
  cellLength = length(cellstring);
  
  for i = 1:nspecies
    
    cellstring{cellLength + i} = strcat('covMat_dot(',num2str(i),') = ',char(covMat_dot_param(i,i)),';');   
  
  end
  
  cellLength = length(cellstring);
  
  counter = 0;
  
  for i = 1:nspecies
  
    for j = i+1:nspecies
      
      counter = counter + 1;
    
      cellstring{cellLength + counter} = strcat('covMat_dot(',num2str(nspecies + counter),') = ',char(covMat_dot_param(i,j)),';');
      
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

    for i = 1:length(xinit_no_sensitivities)

        cellstring{2+i} = strcat('xinit(',num2str(i),') = ',num2str(xinit_no_sensitivities(i)),';');

    end
    
    cellLength = length(cellstring);
    
    functionName = strcat('MRE_',circuitName);
    
    cellstring{cellLength + 1} = 'n = simLength/10;';
    cellstring{cellLength + 2} = 'tspan = linspace(0,simLength,n);';
    cellstring{cellLength + 3} = strcat('[timepoints,Xvals] = ode15s(@',functionName,',tspan,xinit);');
    
    for i = 1:length(cellstring)
        
        temp = cell2mat(cellstring(1,i));
        fprintf(fid,[temp '\n']);
        
    end
    
    fclose(fid);
    
  
  %% Create a file 'MRE_circuitName_.m' to be used by the ODE solver.
  %% Contains S, macroF and covMat_dot. 
  
  str1 = str2mat('touch MRE_unparam_');
  str2 = str2mat(circuitName);
  str3 = str2mat('.m');

  sysCommand = [str1,str2,str3];
    
  system(sysCommand);
  
  filename = strcat('MRE_unparam_',circuitName,'.m');
  
  fid = fopen(filename,'w');
  
  cellstring = {};
  
  cellstring{1} = strcat('function dxdt = MRE_unparam_',circuitName,'(t,x,y)');
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

    str1 = str2mat('touch LNAsim_unparam_');
    str2 = str2mat(circuitName);
    str3 = str2mat('.m');

    sysCommand = [str1,str2,str3];

    system(sysCommand);

    filename = strcat('LNAsim_unparam_',circuitName,'.m');

    fid = fopen(filename,'w');

    cellstring = {};

    cellstring{1} = strcat('function [timepoints,Xvals] = LNAsim_unparam_',circuitName,'(y,simLength)');
    cellstring{2} = '';

    for i = 1:length(xinit_no_sensitivities)

        cellstring{2+i} = strcat('xinit(',num2str(i),') = ',num2str(xinit_no_sensitivities(i)),';');

    end
    
    cellLength = length(cellstring);
    
    cellstring{cellLength + 1} = 'n = simLength/10;';
    cellstring{cellLength + 2} = 'tspan = linspace(0,simLength,n);';
    cellstring{cellLength + 3} = strcat('[timepoints,Xvals] = ode15s(@(t,x) MRE_unparam_',circuitName,'(t,x,y),tspan,xinit);');
    
    for i = 1:length(cellstring)
        
        temp = cell2mat(cellstring(1,i));
        fprintf(fid,[temp '\n']);
        
    end
    
    fclose(fid);
    
    %% Create a file MRE_sensitivities_circuitName.m, used by ODE solver, 
    %% which computes one timestep of the LNA equations and the parameter sensitivities
    
    str1 = str2mat('touch MRE_sensitivities_');
    str2 = str2mat(circuitName);
    str3 = str2mat('.m');
    
    sysCommand = [str1,str2,str3];
    
    system(sysCommand);
    
    filename = strcat('MRE_sensitivities_',circuitName,'.m');
    
    fid = fopen(filename,'w');
    
    cellstring = {};
     
    cellstring{1} = strcat('function dxdt = MRE_sensitivities_',circuitName,'(t,x,y)');
    cellstring{2} = '';
    cellstring{3} = strcat('nspecies = ',num2str(nspecies),';');
    cellstring{4} = strcat('S = zeros(nspecies,',num2str(nrates),');');
    cellstring{5} = '';
    
    cellLength = length(cellstring);
  
      nElementsS = nspecies*nrates;

      % Write S to file

      for k = 1:nElementsS

          cellstring{cellLength+k} = strcat('S(',num2str(k),') = ',num2str(S(k)),';');

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
      
      cellLength = length(cellstring);
      counter = 0;
      for i = 1:dimSensitivities(1)
          
          for j = 1:dimSensitivities(2)
              
              counter = counter + 1;
              cellstring{cellLength + counter} = strcat('Y_dot(',num2str(counter),') = ',char(Y_dot(i,j)),';');
              
          end
          
      end
          
      
      % 'nterms' is the number of terms in the ODE, the first length(x) are the
      % species for the MRE, and the next (length(x)^2 - length(x))/2 + length(x) are the
      % variances and covariances

      cellLength = length(cellstring);
      
      cellstring{cellLength + 1} = '';
      cellstring{cellLength + 2} = 'nterms = nspecies + (nspecies^2 - nspecies)/2 + nspecies;';
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
      cellstring{cellLength + 19} = '';
      cellstring{cellLength + 20} = strcat('nElementsYdot = ',num2str(nElementsYdot),';');
      cellstring{cellLength + 21} = 'for i = 1:nElementsYdot';
      cellstring{cellLength + 22} = '';
      cellstring{cellLength + 23} = 'dxdt(nterms + i) = Y_dot(i);';
      cellstring{cellLength + 24} = '';
      cellstring{cellLength + 25} = 'end';
      
      for i = 1:length(cellstring)

          temp = cell2mat(cellstring(1,i));
          fprintf(fid,[temp '\n']);

      end

      fclose(fid);
      
      %% Create file LNAsim_sensitivities_circuitName.m
      
      str1 = str2mat('touch LNAsim_sensitivities_');
    str2 = str2mat(circuitName);
    str3 = str2mat('.m');

    sysCommand = [str1,str2,str3];

    system(sysCommand);

    filename = strcat('LNAsim_sensitivities_',circuitName,'.m');

    fid = fopen(filename,'w');

    cellstring = {};

    cellstring{1} = '% Function that outputs a given objective function "obj" at steady-state.  For use in optimisation routines.';
    cellstring{2} = '';
    cellstring{3} = strcat('function [obj] = LNAsim_sensitivities_',circuitName,'(y)');
    cellstring{4} = '';
    cellstring{5} = 'simLength = 100;';

    for i = 1:length(xinit_with_sensitivities)

        cellstring{5+i} = strcat('xinit(',num2str(i),') = ',num2str(xinit_with_sensitivities(i)),';');

    end
    
    cellLength = length(cellstring);
    
    functionName = strcat('MRE_sensitivities_',circuitName);
    
    cellstring{cellLength + 1} = 'n = simLength/10;';
    cellstring{cellLength + 2} = 'tspan = linspace(0,simLength,n);';
    cellstring{cellLength + 3} = 'spanLength = length(tspan);';
    cellstring{cellLength + 4} = 'iteration = 0;';
    cellstring{cellLength + 5} = 'ssCondition = 0;';
    cellstring{cellLength + 6} = 'tempXvals = xinit;';
    cellstring{cellLength + 7} = '';
    cellstring{cellLength + 8} = 'while(ssCondition == 0)';
    cellstring{cellLength + 9} = '';
    cellstring{cellLength + 10} = 'iteration = iteration + 1;';
    cellstring{cellLength + 11} = '';
    cellstring{cellLength + 12} = 'if(mod(iteration,50)==0)';
    cellstring{cellLength + 13} = 'iteration';
    cellstring{cellLength + 14} = 'end';
    cellstring{cellLength + 15} = '';
    cellstring{cellLength + 16} = strcat('[timepoints,Xvals] = ode15s(@(t,x) MRE_sensitivities_',circuitName,'(t,x,y),tspan,tempXvals);');
    cellstring{cellLength + 17} = '';
    cellstring{cellLength + 18} = 'newXvals = Xvals(spanLength,:);';
    cellstring{cellLength + 19} = 'absdiff = abs(newXvals-tempXvals);';
    cellstring{cellLength + 20} = '';
    cellstring{cellLength + 21} = 'if(sum(absdiff)<1e-5)';
    cellstring{cellLength + 22} = 'ssCondition = 1;';
    cellstring{cellLength + 23} = 'break';
    cellstring{cellLength + 24} = 'else';
    cellstring{cellLength + 25} = 'tempXvals = newXvals;';
    cellstring{cellLength + 26} = 'continue';
    cellstring{cellLength + 27} = 'end';
    cellstring{cellLength + 28} = 'end';
    cellstring{cellLength + 29} = '';
    cellstring{cellLength + 30} = '% Introduce expression for objective function here';
    cellstring{cellLength + 30} = 'obj = tempXvals(1);';
    cellstring{cellLength + 31} = '% Optionally introduce gradient here';
    
    for i = 1:length(cellstring)
        
        temp = cell2mat(cellstring(1,i));
        fprintf(fid,[temp '\n']);
        
    end
    
    fclose(fid);
  
end
