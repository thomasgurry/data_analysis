% vector x contains [theta1,theta2,theta3] i.e. rotation angles

function error = rotate_rmsd(x, query_coord, template_coord)

ang = x;

Rx = [1 0 0; 0 cos(ang(1)) -sin(ang(1)); 0 sin(ang(1)) cos(ang(1))];
Ry = [cos(ang(2)) 0 sin(ang(2)); 0 1 0; -sin(ang(2)) 0 cos(ang(2))];
Rz = [cos(ang(3)) -sin(ang(3)) 0 ; sin(ang(3)) cos(ang(3)) 0; 0 0 1];
pre_err = (Rz * (Ry * (Rx * query_coord)) - template_coord).^2;
error = sum(sqrt(sum(pre_err)));