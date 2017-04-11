%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%         Generalised function that returns ODE equations for         %
%             species, var/covariances and sensitivities              %
%              according to Linear Noise Approximation                %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Takes as input the name of the circuit as a string, and the name of the
%%% M file with parameters, rates and stoichiometry matrix, also as a string
%%% (e.g. 'input_sge').  Returns all equations dx/dt in symbolic form.

function [allEquations] = getLNAequations(circuitName)
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
  
  % Replace all 'y1y' (parameter sensitivities) with 'x(N+1)x'
  N = length(symX);
  counter = 0;
  
   for i = 1:nElementsYdot
      
      counter = counter + 1;
      oldvar = strcat('y',num2str(counter),'y');
      newvar = strcat('x',num2str(N+counter),'x');
      
      for j = 1:nElementsYdot
          
          oldterm = char(Y_dot(j));
          newterm = strrep(oldterm,oldvar,newvar);
          Y_dot(j) = sym(newterm);
          
      end
      
  end
  
  %% Substitute parameters by 'k1' by 'y(1)' in rates 'macroFnew' and 
  %% by their numeric value in 'macroFparam' (numeric values specified in
  %% the input file).
  

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
  
  
  %% Substitute all parameters in covMat_dot to 'y(1)' type and create a
  %% covMat_dot_param containing the numeric values of the parameters.
  
  covMat_dot_new = covMat_dot;
  %covMat_dot_param = covMat_dot;
   
  for i = 1:length(param)
    
    oldvar = char(param(i));
    newvar1 = strcat('y(',num2str(i),')');
    newvar2 = strcat(num2str(eval(oldvar)));
    
    for q = 1:nElements
      
      oldcov1 = char(covMat_dot_new(q));
      newcov1 = strrep(oldcov1,oldvar,newvar1);
      covMat_dot_new(q) = sym(newcov1);
      
    end
    
  end
  
  covMat_dot = covMat_dot_new;
  
  toc
  
%% Convert parameters 'k1' to 'y(1)' in 'LNAequations'

LNAequations_new = LNAequations;

for i = 1:length(param)
    
    oldvar = param(i);
    str1 = 'y(';
    str2 = ')';
    newvar = strcat(str1,num2str(i),str2);
    
    for j = 1:length(LNAequations)
        
        oldEquation = char(LNAequations_new(j));
        newEquation = strrep(oldEquation,oldvar,newvar);
        LNAequations_new(j) = sym(newEquation);
        
    end
    
end

LNAequations = LNAequations_new;

%% Write all equations into a single symbolic vector to solve dx/dt = 0.

counter = 0;
allEquations = sym([]);

for i = 1:length(LNAequations)
    
    counter = counter + 1;
    allEquations(counter) = LNAequations(i);
    
end

Ydot_dim = size(Y_dot); 

for i = 1:Ydot_dim(1)
    
    for j = 1:Ydot_dim(2)
    
        counter = counter + 1;
        allEquations(counter) = Y_dot(i,j);
        
    end
    
end
  
end