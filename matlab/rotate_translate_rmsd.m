% vector x contains [x,y,z,theta1,theta2,theta3] i.e. translation vector
% and rotation angles

function error = rotate_translate_rmsd(x, coord1, template_coord)

ang = x(4:6);
xyz = x(1:3).';

% First objective function: the translation vector

ncoords = size(coord1,2);

obj1 = 0;

for i = 1:ncoords

    obj1 = ((coord1(:,i) + xyz) - template_coord(:,i)).^2;

end

obj1 = sum(obj1);

% Second objective function: the rotation angles

Rx = [1 0 0; 0 cos(ang(1)) -sin(ang(1)); 0 sin(ang(1)) cos(ang(1))];
Ry = [cos(ang(2)) 0 sin(ang(2)); 0 1 0; -sin(ang(2)) 0 cos(ang(2))];
Rz = [cos(ang(3)) -sin(ang(3)) 0 ; sin(ang(3)) cos(ang(3)) 0; 0 0 1];
pre_err = (Rz * (Ry * (Rx * coord1)) - template_coord).^2;
obj2 = sum(sqrt(sum(pre_err)));

error = 1000*obj1 + obj2;