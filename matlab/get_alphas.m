%% Function that reads ensemble weights from a file of type 
%% 'vbw_output.txt' as output from the VBW program.

function [alphas] = get_weights(filename)

fid = fopen(filename);
lines = fgetl(fid);
lines = fgetl(fid);
lines = fgetl(fid);
alphas = lines(23:length(lines));
alphas = str2num(alphas);
fclose(fid);

