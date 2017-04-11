%%% Function that plots the results of a LNA simulation.
%%% Takes as inputs the index number of the species to be plotted.

% Get indices of timepoints at which the variance is to be displayed
  
function plotLNA(speciesIndex,t,x,varargin)

isColour = nargin - 3;  % Is colour defined

nvars = size(x,2);

nspecies = -3/2 + sqrt(2*nvars + 9/4);


  indexIncr = length(t)/10; 
  varIndices = [indexIncr:indexIncr:length(t)];
  
  xvals = zeros(1,nvars);
  xValsPlot = zeros(10,nvars);
  tPlot = zeros(1,10);
  counter = 0;
  
  for i = varIndices
      
      counter = counter + 1;
      xvals = x(i,:);
      xValsPlot(counter,:) = xvals;
      tPlot(counter) = t(i);
      
  end
  
 
  %% Plot the resulting timeseries
  
  if (isColour == 0)
  
      plot(t,x(:,speciesIndex))
      hold on
      errorbar(tPlot,xValsPlot(:,speciesIndex),sqrt(xValsPlot(:,speciesIndex+nspecies)),'b+')
      hold off
  
  elseif (isColour == 1)
      
      colourStr = varargin{1};
      errorBarStr = strcat(varargin{1},'+');
      
      plot(t,x(:,speciesIndex),colourStr)
      hold on
      errorbar(tPlot,xValsPlot(:,speciesIndex),sqrt(xValsPlot(:,speciesIndex+nspecies)),errorBarStr)
      hold off
  
  else
      
      fprintf('Invalid number of input arguments.')
      
  end
      
  end