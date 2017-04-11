%% Function that reads ensemble weights from a file of type 
%% 'vbw_output.txt' as output from the VBW program.

function [weights] = get_weights(filename)

fid = fopen(filename);
lines = fgetl(fid);
weights = lines(31:length(lines));
weights = str2num(weights);
fclose(fid);

