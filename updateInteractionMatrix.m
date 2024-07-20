% calculate the interaction matrix according to the epigenetic marks
function inter_matrix = updateInteractionMatrix(params,x,Rouse_matrix)

position_left = x(1,1)+1;
position_right = x(1,2)+1;

bd = zeros(params.monomer_num);
bd(position_left,position_right) = position_left-position_right;
bd(position_right,position_left) = position_left-position_right;

A = (bd ~= 0)& (bd ~= -1) ; % attraction
A = A.*(-params.attraction_coef);
R = (bd ~= 0)& (bd == -1) ; % replusion
R = R.*params.repulsion_coef;
T = A + R;
d = diag(true(1,params.monomer_num));
T(d) = -sum(T,2);
inter_matrix = Rouse_matrix + T;
end