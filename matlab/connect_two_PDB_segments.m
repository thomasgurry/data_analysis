%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Script that places the backbone of 4 aSyn helices on the top of a %%%%
%%%%              four-helix bundle by backbone atom RMSD              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: ferritin helices are extracted from the four helix bundle found in
% ferritin where helix 1 has been translated so that its centre of mass is
% at the origin.

%% Define residues from aSyn helix and which residues we want them on the template

clear all
clear alignby_coord1;
clear alignby_coord2;

first_connecting_residue = 187;
second_connecting_residue = 1;


%% PDB's to use

CsgA_pdb = pdbread('csga_linker_CBD.pdb');
linker_pdb = pdbread('cbd_c_term_with_7_his.pdb');

%% Read coordinates of CsgA

tmp = struct2cell(CsgA_pdb.Model.Atom);
natoms = length(tmp);

% Backbone data
isbackbone = squeeze(strcmp(tmp(2,1,:),'CA')) | squeeze(strcmp(tmp(2,1,:),'C')) |squeeze(strcmp(tmp(2,1,:),'O'))| squeeze(strcmp(tmp(2,1,:),'N'));%strcmp compares stings, outputs 1


    bb_coord{1} = squeeze(cell2mat(tmp(8:10,1,isbackbone))); % fetches coordinates of backbone atoms
    bb_atomtype{1} = squeeze(tmp(2,1,isbackbone)); % name of backbone atom
    bb_resname{1} = squeeze(tmp(4,1,isbackbone)); %the names of the residues to which the backbone atoms belong
    bb_resid{1} = squeeze(cell2mat(tmp(6,1,isbackbone)));%the residue id for which the backbone atoms belong to
    bb_atomnum{1} = squeeze(cell2mat(tmp(1,1,isbackbone)));
    %==================================Side Chain data=============================================================
    sc_coord{1} = squeeze(cell2mat(tmp(8:10,1,~isbackbone))); %Gives the coordinates of the side chain atoms in the form of a table, these are the non (~) backbone atoms.
    sc_atomtype{1} = squeeze(tmp(2,1,~isbackbone));%the name of the sidechain atom
    sc_resname{1} = squeeze(tmp(4,1,~isbackbone)); % The SideChain residue names for which the atoms belong to
    sc_resid{1} = squeeze(cell2mat(tmp(6,1,~isbackbone))); % The SideChain residue ID for which the atoms belong to
    sc_atomnum{1} = squeeze(cell2mat(tmp(1,1,~isbackbone)));

    
%% Read coordinates of linker + Mfp3


    tmp = struct2cell(linker_pdb.Model.Atom);
    
    % Backbone data
    isbackbone = squeeze(strcmp(tmp(2,1,:),'CA')) | squeeze(strcmp(tmp(2,1,:),'C')) |squeeze(strcmp(tmp(2,1,:),'O'))| squeeze(strcmp(tmp(2,1,:),'N'));%strcmp compares stings, outputs 1
    %for the next segment
    bb_coord{2} = squeeze(cell2mat(tmp(8:10,1,isbackbone))); % fetches coordinates of backbone atoms
    bb_atomtype{2} = squeeze(tmp(2,1,isbackbone)); % name of backbone atom
    bb_resname{2} = squeeze(tmp(4,1,isbackbone)); %the names of the residues to which the backbone atoms belong
    bb_resid{2} = squeeze(cell2mat(tmp(6,1,isbackbone)));%the residue id for which the backbone atoms belong to
    bb_atomnum{2} = squeeze(cell2mat(tmp(1,1,isbackbone)));
    %==================================Side Chain data=============================================================
    sc_coord{2} = squeeze(cell2mat(tmp(8:10,1,~isbackbone))); %Gives the coordinates of the side chain atoms in the form of a table, these are the non (~) backbone atoms.
    sc_atomtype{2} = squeeze(tmp(2,1,~isbackbone));%the name of the sidechain atom
    sc_resname{2} = squeeze(tmp(4,1,~isbackbone)); % The SideChain residue names for which the atoms belong to
    sc_resid{2} = squeeze(cell2mat(tmp(6,1,~isbackbone))); % The SideChain residue ID for which the atoms belong to
    sc_atomnum{2} = squeeze(cell2mat(tmp(1,1,~isbackbone)));

%% Compute the correct translation for linker onto CsgA

alignby = {};   
cm = {};

    resid = bb_resid{2};
    alignby{1} = resid>= min(second_connecting_residue) & resid<=max(second_connecting_residue);

    % Calculate the centre of mass of the linker

    cm{1} = mean(bb_coord{2}(:, alignby{1}),2);

    
% Translate and rotate 20 times

    for k = 1:20

        % Calculate the centre of mass of the linker pre-translation
        CsgA_resid = bb_resid{1};
        alignby_CsgA = CsgA_resid>= min(first_connecting_residue) & CsgA_resid<= max(first_connecting_residue);
        cmCsgA = mean(bb_coord{1}(:, alignby_CsgA), 2);
        
        
        % Align centre of masses of aSyn to the respective helices
        linker_translation = cm{1} - cmCsgA;
        
        %% Compute the correct rotation for each helix by minimising RMSD
        
        % Translate all backbone and sidechain atoms for aSyn
        % monomer
        
        for j = 1:size(bb_coord{1},2)
            
            bb_coord{1}(:,j) = bb_coord{1}(:,j) + linker_translation;
            
        end
        
        for j = 1:size(sc_coord{1},2)
            
            sc_coord{1}(:,j) = sc_coord{1}(:,j) + linker_translation;
            
        end
        
        
        % Rotate all backbone and sidechain atoms for each aSyn
        % monomer
        
        
        alignby_coord2 = bb_coord{2}(:, alignby{1});
        alignby_coord1 = bb_coord{1}(:, alignby_CsgA); % contains backbone atoms of the aSyn helix
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    HERE IS THE MODIFICATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make sure they are the same dimension
        alignby_coord2 = alignby_coord2(:,1:3);
                
        angles = fminsearch(@(x)rotate_rmsd(x,alignby_coord1,alignby_coord2), [0.01, 0.01, 0.01]);
        
        ang = angles;
        
        Rx = [1 0 0; 0 cos(ang(1)) -sin(ang(1)); 0 sin(ang(1)) cos(ang(1))];
        Ry = [cos(ang(2)) 0 sin(ang(2)); 0 1 0; -sin(ang(2)) 0 cos(ang(2))];
        Rz = [cos(ang(3)) -sin(ang(3)) 0 ; sin(ang(3)) cos(ang(3)) 0; 0 0 1];
        
        % Rotate helices
        bb_coord{1} = Rz * (Ry * (Rx * bb_coord{1}));
        sc_coord{1} = Rz * (Ry * (Rx * sc_coord{1}));

        k
        
    end


    %% Combine new structure information for CsgA

    all_resid{1} = [bb_resid{1}; sc_resid{1}];
    all_resname{1} = [bb_resname{1}; sc_resname{1}];
    all_coord{1} = [bb_coord{1} sc_coord{1}];
    all_atomtype{1} = [bb_atomtype{1}; sc_atomtype{1}];
    all_atomnum{1} = [bb_atomnum{1}; sc_atomnum{1}];
    
    [atomnum indices] = sort(all_atomnum{1});
    all_atomnum{1} = atomnum(indices);
    all_resid{1} = all_resid{1}(indices);
    all_resname{1} = all_resname{1}(indices);
    all_coord{1} = all_coord{1}(:,indices);
    all_atomtype{1} = all_atomtype{1}(indices);

    %% Combine new structure information for Mfp3

    all_resid{2} = [bb_resid{2}; sc_resid{2}];
    all_resname{2} = [bb_resname{2}; sc_resname{2}];
    all_coord{2} = [bb_coord{2} sc_coord{2}];
    all_atomtype{2} = [bb_atomtype{2}; sc_atomtype{2}];
    all_atomnum{2} = [bb_atomnum{2}; sc_atomnum{2}];
    
    [atomnum indices] = sort(all_atomnum{2});
    all_atomnum{2} = atomnum(indices);
    all_resid{2} = all_resid{2}(indices);
    all_resname{2} = all_resname{2}(indices);
    all_coord{2} = all_coord{2}(:,indices);
    all_atomtype{2} = all_atomtype{2}(indices);
    
    %% Concatenate into new full structure
    
    combined_resid = all_resid{1};
    combined_resname = all_resname{1};
    combined_coord = all_coord{1};
    combined_atomtype = all_atomtype{1};
    combined_atomnum = all_atomnum{1};
    
    isnot_overlap = all_resid{2}~=1;
    combined_resid = [combined_resid; (all_resid{2}(isnot_overlap)+first_connecting_residue-1)];
    combined_resname = [combined_resname; all_resname{2}(isnot_overlap)];
    combined_coord = [combined_coord, all_coord{2}(:,isnot_overlap)];
    combined_atomtype = [combined_atomtype; all_atomtype{2}(isnot_overlap)];
    combined_atomnum = [combined_atomnum; all_atomnum{2}(isnot_overlap)];
    
    
%% Convert HIS to HSE

for i = 1:length(combined_resname)
    
    if(strcmp(combined_resname(i),'HIS'))
        
        combined_resname{i} = 'HSD';
        
    end
    
end
    
%% Write out combined structure PDB for fusion construct

fid = fopen(sprintf('csga_linker_CBD_7his.pdb'),'w+');
fprintf(fid,'REMARK CsgA + linker + CBD + His tag structure\n');

for atomnumber = 1:size(combined_resid,1)
    
    fprintf(fid,'ATOM %6d %4s %3s %5d %10.3f %8.3f %8.3f  1.00  0.00      A\n', atomnumber, char(combined_atomtype(atomnumber)), char(combined_resname(atomnumber)), combined_resid(atomnumber), str2num(num2str(combined_coord(1,atomnumber))), str2num(num2str(combined_coord(2,atomnumber))), str2num(num2str(combined_coord(3,atomnumber))));
    
end

fclose(fid);

