%% Function that reads ensemble weights from a file of type 
%% 'ensemble_summary.txt' as output from the VBW program.

function [weights,alphas,files] = parse_ensemble_summary(filename)

fid = fopen(filename);
file_contents = textscan(fid,'%s %s %s','HeaderLines',1);
fclose(fid)

% Figure out how many PDB structures there are
npdbs = find(strcmp(file_contents{1},'Scaling')) - 1;

files = file_contents{1}(1:npdbs);
weights = file_contents{2}(1:npdbs);
alphas = file_contents{3}(1:npdbs);

weights = str2num(char(weights));
alphas = str2num(char(alphas));