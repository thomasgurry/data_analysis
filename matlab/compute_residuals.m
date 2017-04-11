%% Function that computes and plots the residuals between experimental 
%% and predicted data.  Takes as input the name of the ensemble_data.dat file.

function [residuals_matrix] = compute_residuals(datafile)

fid = fopen(datafile);
l = fgetl(fid);
file_contents = textscan(fid,'%d %s %s %f %f');

% Get logical vectors for each measurement type

meas_types = {'HA','H','N','CA','CB','HN-RDC','CCA-RDC','CN-RDC'};
labels = file_contents{3};
residue_numbers = file_contents{1};
nresidues = max(residue_numbers);

matrix_of_logicals = zeros(length(labels),length(meas_types));

for i = 1:length(labels)
    
    for j = 1:length(meas_types)
        
        matrix_of_logicals(i,j) = strcmp(labels{i},meas_types(j));
        
    end
    
end

% Compute residuals for each measurement type and store 
% in a matrix to be output.  Note that a non-existent residual
% i.e. one with no experimental data for that residue, is computed 
% as NA

residuals_matrix = NaN(nresidues,length(meas_types));
exp_matrix = NaN(nresidues,length(meas_types));
pred_matrix = NaN(nresidues,length(meas_types));
exp_data = file_contents{4};
pred_data = file_contents{5};

for i = 1:length(meas_types)

    logical_vect = logical(matrix_of_logicals(:,i));
    resIDs = residue_numbers(logical_vect);
    residuals = exp_data(logical_vect) - pred_data(logical_vect);
    residuals_matrix(resIDs,i) = residuals;
    exp_matrix(resIDs,i) = exp_data(logical_vect);
    pred_matrix(resIDs,i) = pred_data(logical_vect);
    
end

% Plot residuals for each measurement type on separate plots.
    
total_errors = [0.23, 0.49, 2.43, 0.98, 1.1, 1.16, 0.4, 0.5, 0.3]; % total^2 = (exp^2 + pre^2)

% Plot residuals
figure
for i = 1:(length(meas_types))
    
    subplot(2,4,i)
    plot([1:nresidues],residuals_matrix(:,i),'*')
    xlabel('Residue number')
    ylabel('Experimental - Predicted')
    title(strcat('Measurement type: ',meas_types{i}))
    zeroline = line([1,nresidues],[0,0])
    set(zeroline,'Color','k')
    errorline1 = line([1,nresidues],[total_errors(i),total_errors(i)]);
    errorline2 = line([1,nresidues],[-total_errors(i),-total_errors(i)]);
    set(errorline1,'Color','r')
    set(errorline2,'Color','r')
    
end

% Plot both experimental and predicted data on the same plots

figure
for i = 1:(length(meas_types))
    
    subplot(2,4,i)
    plot([1:nresidues],exp_matrix(:,i),'-b')
    hold on
    plot([1:nresidues],pred_matrix(:,i),'-r')
    xlabel('Residue number')
    ylabel('Measurement value')
    title(strcat('Measurement type: ',meas_types{i}))
    zeroline = line([1,nresidues],[0,0])
    set(zeroline,'Color','k')
    legend('Experimental','Predicted')
    hold off
    
end

% Plot experimental on one axis and predicted on another

maxCalpha = max(pred_matrix(:,4))+5;
minCalpha = min(pred_matrix(:,4))-5;
maxRDC = max(pred_matrix(:,6))+5;
minRDC = min(pred_matrix(:,6))-5;

figure
subplot(1,2,1)
plot(exp_matrix(:,4),pred_matrix(:,4),'.')
hold on
xlabel('Experimental C\alpha chemical shift (ppm)')
ylabel('Ensemble C\alpha chemical shift (ppm)')
axis([minCalpha maxCalpha minCalpha maxCalpha])
line([minCalpha,maxCalpha],[minCalpha,maxCalpha])
hold off

subplot(1,2,2)
plot(exp_matrix(:,6),pred_matrix(:,6),'.')
hold on
xlabel('Experimental HN-RDC (Hz)')
ylabel('Ensemble HN-RDC (Hz)')
axis([minRDC maxRDC minRDC maxRDC])
line([minRDC,maxRDC],[minRDC,maxRDC])
hold off

end