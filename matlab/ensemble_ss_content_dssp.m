% Computes secondary structure content of an ensemble.  Must be called from
% a folder containing dssp output text files and no other text files.

% Input weights, alphas and files from parse_ensemble_summary()

function [helix_assignments,strand_assignments] = ensemble_ss_content_dssp(weights,alphas,summary_files,varargin)

if(size(varargin,2)>0)
    isplot = 0;
else
    isplot = 1;
end


nresidues = 140;

%% Read in dssp ID file

% Check if it exists!!
files = dir('*.txt');
nfiles = length(files);
filenames = {};
cond = 0;

for i = 1:nfiles
    if(strcmp(files(i).name,'dssp_IDs.txt') == 0)
        filenames{i} = files(i).name;
    else
        cond = 1;
    end
end

if(cond == 0)
    error('No dssp_IDs.txt file found!!\n')
else
    nfiles = nfiles - 1;
end

% Read in filenames in order
dssp_IDs_fid = fopen('dssp_IDs.txt','r');
dssp_IDs_contents = textscan(dssp_IDs_fid,'%d %s');
dssp_indices = dssp_IDs_contents{1};
dssp_filenames = {};
fclose(dssp_IDs_fid);

for i = 1:length(dssp_indices)
    
    filename = dssp_IDs_contents{2}{i}(38:length(dssp_IDs_contents{2}{i}));
    dssp_filenames{i} = filename(1:(length(filename)-4));
    
end

%% Find dssp indices of filenames in the order of the summary_files (same
%% order as the weights).

dssp_inds_in_order = [];

% Check that summary_files (from ensemble_summary.txt file) and
% dssp_filenames have the same length.

if(length(summary_files) ~= length(dssp_filenames))
    error('Different number of dssp files and weights!!\n')
end

for i = 1:length(summary_files)
    % Find dssp ID summary_file(i)
    ind = strcmp(dssp_filenames,summary_files(i));
    dssp_inds_in_order(i) = dssp_indices(ind);
end

%% Initialise
helix_assignments = zeros(140,nfiles);
strand_assignments = zeros(140,nfiles);

%% Loop through files and obtain secondary structure assignments
for i = 1:length(dssp_inds_in_order)
    
    filename = strcat('dssp',num2str(dssp_inds_in_order(i)),'.txt');
    %fid = fopen(files(i).name);
    i
    [helix_content, strand_content, nres] = dssp_parser(filename);
    
    helices = zeros(560,1);
    strands = zeros(560,1);
    helices(1:length(helix_content)) = helix_content;
    strands(1:length(strand_content)) = strand_content;
    
    % Check number of alpha-synuclein chains involved and if more than one,
    % average content across the molecule
    
    nchains = nres/140;
    
    if(nchains > 1)
        helices = reshape(helices,140,4);
        strands = reshape(strands,140,4);
        helices = mean(helices(:,1:nchains),2);
        strands = mean(strands(:,1:nchains),2);
        helix_assignments(:,i) = helices;
        strand_assignments(:,i) = strands;
    else
        helix_assignments(:,i) = helices(1:140);
        strand_assignments(:,i) = strands(1:140);
    end
        
end


%% Analyse ensemble content

if(weights == 1)
    weights = (1/nfiles)*ones(nfiles,1);
end

if(isplot == 1)
    
    % Plot distribution of helical content
    figure
    hist(sum(helix_assignments)/140)
    xlabel('Fractional helical content')
    ylabel('Frequency')
    title('Structural library helical content')
    
    % Plot distribution of strand content
    figure
    hist(sum(strand_assignments)/140)
    xlabel('Fractional strand content')
    ylabel('Frequency')
    title('Structural library strand content')

end

%% Plot ensemble helical and strand contents

helical_content = zeros(140,1);
strand_content = zeros(140,1);

if(size(weights,1) == 1)
    weights = weights.'; % make weights vector a column vector
end

for i = 1:140
    
    residue_helical_contents = helix_assignments(i,:);
    residue_strand_contents = strand_assignments(i,:);
    helical_content(i) = helix_assignments(i,:)*weights;
    strand_content(i) = strand_assignments(i,:)*weights;
    
end

% Calculate 95% confidence intervals

alpha0 = sum(alphas);

helical_vars = zeros(140,1);
strand_vars = zeros(140,1);

for n = 1:140
    
    helicalvar = 0;
    strandvar = 0;
    
    for i = 1:140
        
        for j = 1:140
            
            if(i == j)
                kroneckerDelta = 1;
            else
                kroneckerDelta = 0;
            end
            
            helicalvar = helicalvar + helix_assignments(n,i)*helix_assignments(n,j)*(alphas(i)*alpha0*kroneckerDelta - alphas(i)*alphas(j));
            strandvar = strandvar + strand_assignments(n,i)*strand_assignments(n,j)*(alphas(i)*alpha0*kroneckerDelta - alphas(i)*alphas(j));

        end
 
    end
    
    helical_vars(n) = (helicalvar/((alpha0^2)*(alpha0 + 1)));
    strand_vars(n) = (strandvar/((alpha0^2)*(alpha0 + 1)));
    
end

helical_CIs = sqrt(helical_vars)*1.96*1.54;
strand_CIs = sqrt(strand_vars)*1.96*1.54;

if(isplot == 1)

    figure
    plot([1:140],helical_content,'-r')
    hold on
    plot([1:140],strand_content,'-b')
    xlabel('Residue number')
    ylabel('Fractional content')
    title('Ensemble secondary structure propensities')
    legend('Helix','Strand')
    
    figure
    errorbar([1:140],helical_content,helical_CIs,'-r')
    hold on
    errorbar([1:140],strand_content,helical_CIs,'-b')
    xlabel('Residue number')
    ylabel('Fractional content')
    title('Ensemble secondary structure propensities')
    legend('Helix','Strand')
    
    figure
    plot([1:140],helical_content+helical_CIs,'-r')
    hold on
    plot([1:140],helical_content-helical_CIs,'-r')
    plot([1:140],strand_content+strand_CIs,'-b')
    plot([1:140],strand_content-strand_CIs,'-b')
    xlabel('Residue number')
    ylabel('Fractional content')
    title('Ensemble library secondary structure propensities')
    legend('Helix','Helix','Strand','Strand')

end

% Plot unweighted helical and strand contents

unweighted_helical_content = zeros(140,1);
unweighted_strand_content = zeros(140,1);

for i = 1:140
    
    residue_helical_contents = helix_assignments(i,:);
    residue_strand_contents = strand_assignments(i,:);
    unweighted_helical_content(i) = mean(helix_assignments(i,:));
    unweighted_strand_content(i) = mean(strand_assignments(i,:));
    
end

if (isplot == 1)

    figure
    plot([1:140],unweighted_helical_content,'-r')
    hold on
    plot([1:140],unweighted_strand_content,'-b')
    xlabel('Residue number')
    ylabel('Fractional content')
    title('Unweighted library secondary structure propensities')
    legend('Helix','Strand')

end
