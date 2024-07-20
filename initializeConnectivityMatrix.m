function R = initializeConnectivityMatrix(params)
R                                        = zeros(params.monomer_num);
R(diag(true(1,params.monomer_num-1),1))  = -1;  % super diagonal
R(diag(true(1,params.monomer_num-1),-1)) = -1;  % sub-diagonal

% Multiply by the spring constant values for heterogeneous polymer
R    = R.*params.spring_const;
% Sum rows to get the diagonal element value
d    = diag(true(1,params.monomer_num));
R(d) = -sum(R,2);
end