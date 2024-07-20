function params = parametersChromatinDyn(monomer_num)
    params.monomer_num    = monomer_num ; % nucleosome number
    params.dt                = 0.01;  % time step
    params.dimension         = 3; % spatial dimension
    params.friction_coef     = 1; % friction coefficient [Newton][sec]/[mu m]
    params.diffusion_const   = 1; % diffusion constant [mu m]^2 /[sec]
    params.b                 = sqrt(3); % STD of distance between adjacent monomers
    params.spring_const      = (params.dimension.*params.diffusion_const./params.b^2) .*ones(params.monomer_num); % spring constant
    params.attraction_coef   = params.spring_const/1000; % attraction interaction constant
    params.repulsion_coef    = 0; % repulsion interaction constant
    params.factor            = sqrt(2*params.diffusion_const*params.dt);
end
