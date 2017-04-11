% vector x contains [x,y,z,theta1,theta2,theta3] i.e. translation vector
% and rotation angles

function error = rotate_rmsd_angles(x, coord1, template_coord)

ang = x;

% Second objective function: the rotation angles

Rx = [1 0 0; 0 cos(ang(1)) -sin(ang(1)); 0 sin(ang(1)) cos(ang(1))];
Ry = [cos(ang(2)) 0 sin(ang(2)); 0 1 0; -sin(ang(2)) 0 cos(ang(2))];
Rz = [cos(ang(3)) -sin(ang(3)) 0 ; sin(ang(3)) cos(ang(3)) 0; 0 0 1];
pre_err = (Rz * (Ry * (Rx * coord1)) - template_coord).^2;
error = sum(sqrt(sum(pre_err)));