% vector x contains [x,y,z,theta1,theta2,theta3] i.e. translation vector
% and rotation angles

function obj1 = translate_rmsd(x, coord1, template_coord)

xyz = x(1:3).';

% First objective function: the translation vector

ncoords = size(coord1,2);

obj1 = 0;

for i = 1:ncoords

    obj1 = ((coord1(:,i) + xyz) - template_coord(:,i)).^2;

end

obj1 = sum(obj1);