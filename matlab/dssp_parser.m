function [ helix_assignments, strand_assignments, nresidues ] = dssp_parser( file )

% Parse DSSP output file into two binary strings for helical and strand
% content

helix_assignments = zeros(560,1);
strand_assignments = zeros(560,1);

fid = fopen(file);
dssp_res = textscan(fid,'%d %s');
resnums = dssp_res{1};
maxres = max(resnums);
if(maxres < 141)
    nresidues = 140;
elseif(maxres > 140 & maxres < 281)
    nresidues = 280;
elseif(maxres > 280 & maxres < 421)
    nresidues = 420;
elseif(maxres > 420)
    nresidues = 560;
end
    
assignments = dssp_res{2};

for i = 1:length(assignments)

    if(strcmp(assignments(i),'H') | strcmp(assignments(i),'G') | strcmp(assignments(i),'I'))
        helix_assignments(resnums(i)) = 1;
    elseif(strcmp(assignments(i),'B') | strcmp(assignments(i),'E') | strcmp(assignments(i),'I'))
        strand_assignments(resnums(i)) = 1;
    end
    
end

fclose(fid);

end

