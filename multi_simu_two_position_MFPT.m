% % Simulate the one-dimensional case of the mean first arrival time (MFPT) required for cohesin subunits to diffuse to the CTCF boundary
% % physical barrier 
clear;clc

% % Parameter settings
N      = 100;    % Lattice Number
k_diff = 35;     % the subunit diffuses rate 
d      = 0.05;   
d_off  = d/2;    % the subunit dissociation rate
%  epsilon  = 1.15;   %  kphys = k_diff*exp(-epsilon)
epsilon = 0.0001:0.75:6.0001;

right_loading_site = 25;  % right subunit loading position
% barrier_site       = 7;  % the physical barrier position
barrier_site = [3    12    21    30   39    48       57    66    75    84  93    99];

simu_MFPT=zeros(length(barrier_site),length(epsilon)); % pre-allocating matrix

for i =1:length(barrier_site)
    for j = 1:length(epsilon)
        simu_MFPT(i,j) = a_simu_two_position_MFPT(N,k_diff,d_off,right_loading_site,barrier_site(i),epsilon(j));
        disp(i)
    end
end

save('simu_N=100LS=(24,25)k=35d=0.05beta=0.05.mat',"simu_MFPT","d",'k_diff',"d_off",'barrier_site','epsilon')

%%  Multiple simulations to calculate  MFPT
function avg_t_absorbed =a_simu_two_position_MFPT(N,k_diff,d_off,right_loading_site,barrier_site,epsilon)
tic
num_runs         = 1e5;                 % the number of simulation
t_absorbed_array = zeros(num_runs, 1);  % create an empty array to save each run of t_absorbed

parfor j = 1:num_runs
%     3D simulation 
    params               = parametersChromatinDyn(N+1);
   [~,~,~, t_absorbed,~] = Gillespie_two_position_3D(params,N,k_diff,d_off,right_loading_site,barrier_site,epsilon);
   t_absorbed_array(j)   = t_absorbed;
   disp(j)
end
toc

% calculate the average value of t_absorbed
New_t_absorbed_array = t_absorbed_array(t_absorbed_array ~= 0);
avg_t_absorbed       = mean(New_t_absorbed_array); 
end
