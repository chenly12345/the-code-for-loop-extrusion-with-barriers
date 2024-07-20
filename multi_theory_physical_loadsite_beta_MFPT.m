% Theoretical data for different epsilon values and physical barrier positions
clear;clc

num_N    = 100;  % Length of lattice numbering
params.N = num_N - 1; % Lattice numbering from 0 to N, matrix size N*(N+1)/2
params.k = 35; % Cohesin transition rate
params.d = 0.05;  % Cohesin dissociation rate (k_off in text)
TT       = 1200; % Integration interval for ODE solving

% right_load_site =2:2:50;
right_load_site = 25; % Right loading site for cohesin

% indices_b       = 72; 
indices_b = 1:num_N -1; % Barrier positions excluding right loading site
indices_b([(right_load_site - 1), right_load_site]) = [];  % Removing obstacle positions

num_site = ((right_load_site + 1) * right_load_site) / 2; % Initial cohesin loading positions

epsilon = 0.001:0.075:6; % Range of beta values
% epsilon = 1.15;

T_MFPT_barrier=zeros(length(epsilon),length(indices_b));

% Loop over epsilon values and barrier positions
for i =1:length(epsilon)
    for j = 1:length(indices_b)
        % If the obstacle site is at the right end of the initial position, P0 does not change. 
        % If the obstacle site is at the left end of the initial position, P0 moves forward two positions.
        if right_load_site <indices_b(j)
           n= (right_load_site+1)*right_load_site/2;
        else 
          n= (right_load_site-1)*(right_load_site)/2;
        end
 
        % Calculating mean first passage time with physical barrier
        T_MFPT_barrier(i,j) =single_theory_physical_barrier(params.N,k,params.d,epsilon(i),n,indices_b(j),TT);
    end 
    disp(i)
end

% Calculating theoretical mean first passage time without barrier
T_theory_Nobarrier =MFPT_No_barrier(num_N,k,params.d,num_site,TT);
% Saving results to file
save('N=100LS=(24,25)k=35d=0.05barrier_site=2_99beta=0_6.mat','T_MFPT_barrier','T_theory_Nobarrier','indices_b',"epsilon",'right_load_site')


% Generating and plotting heat map of relative differences
loading_value = nan(length(epsilon),2);
T_value       = (T_MFPT_barrier-T_theory_Nobarrier)'./T_theory_Nobarrier;

% Trimming data for large heatmap import into CDR
epsilon =epsilon(:,1:2:end-3);T_value_last=T_value_last(1:2:end-3,:);
indices_b1    = 1:99;
pcolor(indices_b1./100,epsilon,T_value_last)
% pcolor(indices_b./100,beta,(T_MFPT_barrier-T_theory_Nobarrier)'./T_theory_Nobarrier)
shading interp
xlabel('Barrier site/N')
ylabel('\epsilon')
set(gcf,'Position',[100,100,400,300])

figure
A=(T_MFPT_barrier-T_theory_Nobarrier)./T_theory_Nobarrier;
a=15;
plot(indices_b(1:right_load_site-2)./100,A(1:right_load_site-2,a),"LineWidth",1,"Color",'#0072BD') % Blue
hold on 
plot(indices_b(right_load_site-1:97)./100,A((right_load_site-1):97,a),"LineWidth",1,"Color",'#0072BD')

a=27;
plot(indices_b(1:right_load_site-2)./100,A(1:right_load_site-2,a),"LineWidth",1,"Color",'#FF0000') % Red
plot(indices_b(right_load_site-1:97)./100,A(right_load_site-1:97,a),"LineWidth",1,"Color",'#FF0000')

a=45;
plot(indices_b(1:right_load_site-2)./100,A(1:right_load_site-2,a),"LineWidth",1,"Color",'#77AC30') %Green
plot(indices_b(right_load_site-1:97)./100,A(right_load_site-1:97,a),"LineWidth",1,"Color",'#77AC30')

plot([0 1],[0 0],'--',"LineWidth",1,"Color",'#7E2F8E') %Purple

xlabel('Barrier site/N')
ylabel('\Delta\tau')
set(gcf,'Position',[100,100,400,300])
legend('\epsilon=0.05','','\epsilon=1.15','','\epsilon=2.25','','')
title('Loadingsite 0.25')