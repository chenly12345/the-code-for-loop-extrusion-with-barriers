function T_MFPT = single_theory_physical_barrier(a,b,c,epsilon,n,ind_b1,TT)
    params.N = a; %% Lattice number from 0 to N, matrix size is N*(N+1)/2 -1 
    params.k = b;  %% k is the cohesin transition rate
    params.d = c;  %% cohesin dissociation rate, referred to as k_off in the text
    P_0      = zeros(params.N*(params.N+1)/2-1,1);%% Initial condition 
    P_0(n)   = 1; % Initial value
    delta    = 0.1;
    T        = 0:delta:TT;
   indices_b = ind_b1 ; % Barrier site 

 %% Generate transfer rate matrix
    % Step 1
    % Rate matrix for two cohesin monomers sliding on a chromatin lattice without obstacles, no dissociation rate
    A      = zeros(sum(1:params.N)+params.N);% Initialize a matrix of size N*(N+1)/2 + N
    A(1,2) = params.k; % Lattice boundary
    temp   = 0;%% Initial value
    
    for i = 2:params.N
        % Middle part
        B = A(2+temp:2+temp+i-1,2+temp:2+temp+i-1);%% Extract submatrix B from A, size i*i
        B(diag(true(1,size(B,1)-1),1)) = params.k;  %% Set values above the main diagonal to k
        A(2+temp:2+temp+i-1,2+temp:2+temp+i-1) = B;%% Put B back into A
        
        % Above the main diagonal
        C = A(2+temp:2+temp+i-1,2+temp+i:2+temp+2*i-1);%% Extract submatrix C of the same size as B
        C(diag(true(1,size(C,1)))) =params.k;%% Set values on the diagonal to k
        A(2+temp:2+temp+i-1,2+temp+i:2+temp+2*i-1) = C;
        temp = temp + i; 
    end
    
    A      = A + A'; % Fill the lower triangle of the main diagonal
    A_d    = diag(true(1, size(A, 1))); % Main diagonal of matrix A
    A(A_d) = -sum(A, 2); % Set elements of the main diagonal to the negative sum of each row
    T_v    = A(1:sum(1:params.N), 1:sum(1:params.N)); %% Transition rate matrix for N lattices
    
    % Step 2 
    % Considering absorbing boundary:
    % Lattice number (0, ) decreases by one k; (1, ) decreases one inflow;
    index    = 1: params.N-1; % Specify positions to be replaced
    index1   = cumsum(index) + 1;% Corresponding number for (0, )
    index1(1)= 2;
    
    for i = 1:length(index1)-1
        T_v(index1(i),index1(i)) = -2*params.k;
        T_v(index1(i)+1,index1(i)) = 0;
    end
     
    T1_diag         = diag(eye(params.N));
    T1_diag(1)      = 3*params.k; % (0, N) outflow probability is 0
    T1_diag(end)    = params.k;   % (N-1, N) outflow probability is -k
    T1_diag(2:end-1)= 2*params.k; % ( , N) the rest outflow probability is -2k
    T1              = diag(T1_diag);
    T1(2,1)         = -params.k;
    
    % Expand T1 matrix to (2N-1) dimensions
    T2         = zeros((2*params.N)-1);
    T2(params.N:(2*params.N)-1,params.N:(2*params.N)-1) = T1;
    indices    = [1:params.N-1;params.N:2*params.N-2]';
    new_values = ones(params.N-1,1)*-params.k;
    
    for i = 1:size(indices, 1)
        T2(indices(i,1), indices(i,2)) = new_values(i);
    end
    
    % Pad T2 matrix to the same size as T_v
    T3  = padarray(T2, [sum(1:params.N)-2*params.N+1; sum(1:params.N)-2*params.N+1], 'pre');
    T_k = T_v+T3;

    % Step 3 there are physical barrier ;  the rate k*exp(-beta)
    ind_b    = indices_b-1;
    b        = 0:1:ind_b-1;
    site_by  = sort(ind_b*(ind_b+1)/2-b) ; 
    site_bx  = site_by+ind_b ; 
 
    T_k(sub2ind(size(T_k), site_bx,site_by)) = params.k*exp(-epsilon);
    
    if indices_b ~= params.N
        T_k(sub2ind(size(T_k), site_by,site_bx)) = params.k*exp(-epsilon);
    end
    
    bb = ind_b+1:params.N-1;
    
    if  isempty(site_bx)
        site_b2y=1 +cumsum(bb); 
    else
        site_b2y=site_bx(end)+1 +cumsum(bb); 
    end
    
    site_b2x = site_b2y+1;
    
    if indices_b ~= 1
        T_k(sub2ind(size(T_k), site_b2x,site_b2y)) = params.k*exp(-epsilon);
    end
    
    T_k(sub2ind(size(T_k), site_b2y,site_b2x)) = params.k*exp(-epsilon);
    
    T_k     = T_k-diag(diag(T_k));
    diag_Tk = sum(T_k);
    T_k     = T_k-diag(diag_Tk);
    ind_b   = indices_b-1;
    b       = 0:1:ind_b-1;
    site_by = sort(ind_b*(ind_b+1)/2-b) ; 
    site_bx = site_by+ind_b ; 
    
    T_k(sub2ind(size(T_k), site_bx,site_by)) = params.k*exp(-epsilon);
    
    if indices_b ~= params.N
        T_k(sub2ind(size(T_k), site_by,site_bx)) = params.k*exp(-epsilon);
    end
    
    bb = ind_b+1:params.N-1;
    
    if  isempty(site_bx)
        site_b2y=1 +cumsum(bb); 
    else
        site_b2y=site_bx(end)+1 +cumsum(bb); 
    end
    
    site_b2x = site_b2y+1;
    
    if indices_b ~= 1
        T_k(sub2ind(size(T_k), site_b2x,site_b2y)) = params.k*exp(-epsilon);
    end
    
    T_k(sub2ind(size(T_k), site_b2y,site_b2x)) = params.k*exp(-epsilon);
    
    T_k     = T_k-diag(diag(T_k));
    diag_Tk = sum(T_k);
    T_k     = T_k-diag(diag_Tk);
 
  % Step 4
    % Considering cohesin dissociation on the lattice, dissociation rate is set to d;
    I   = eye(size(T_k));%% Generate an identity matrix of the same size as the matrix without dissociation
    T_d = params.d.*I;   %% Dissociation rate is d 
    
    % If one subunit is at the boundary, the dissociation rate is d/2.
    NN       = params.N;
    d_sigle1 = zeros(1,NN);
    
     for i = 1:NN
         d_sigle1(1,i) = i*(i-1)/2+1;
     end
    
    d_sigle2 = d_sigle1(end):NN*(NN+1)/2;
    d_sigle  = [d_sigle1(1:end),d_sigle2(2:end)];
    
    for i = 1:length(d_sigle)
        T_d(d_sigle(i),d_sigle(i)) =params.d/2;
    end
    
    glob_T=T_k-T_d; %% Transition rate matrix for cohesin diffusion and dissociation on the lattice
    
    % Step 4
    % Remove boundary (0, N), and renumber
    glob_TT = glob_T([1:params.N*(params.N-1)/2 params.N*(params.N-1)/2+2:end],[1:params.N*(params.N-1)/2 params.N*(params.N-1)/2+2:end]);
    glob_T1 = sparse(glob_TT);
    
    %% Solve main equation
    [~,Pt] = solve_ode(glob_T1,P_0,delta,TT);
    P_t    = Pt';
    % P1     = sum(P_t); % Matrix P_t, rows are time t, columns are probability P
    
%% Calculate Mean First Passage Time (MFPT)
    if params.d==0
        N_t        = sum(P_t) ;% Survival probability without dissociation
        T_MFPT     = trapz(T, N_t);
    else

    d_vector=params.d*ones(params.N*(params.N+1)/2,1);
    d_vector(d_sigle(1:end)) = params.d/2;
    d_vector(params.N*(params.N-1)/2+1) = [];

      Fun_int = d_vector.*P_t; % Integrand function
      P_int   = trapz(T, Fun_int,2);
      P_off   = sum(P_int);

    if indices_b == params.N
        a_s=[params.k*exp(-epsilon);params.k]; 
        elseif indices_b  == 1
        a_s=[params.k;params.k*exp(-epsilon)];
        else
        a_s=params.k*ones(2,1); 
    end

        % Probability flow at the boundary: state transition rate * probability at the boundary
        J_s         = a_s.*P_t([(params.N-1)*(params.N-2)/2+1,params.N*(params.N-1)/2+1],:);%(0，N-1)and(1,N)
        Fun_sJ      = T.*J_s;%被积函数
        P_sJ_result = trapz(T, Fun_sJ,2);
        P_total_sJ  = sum(P_sJ_result);
        
        % Calculate the first passage time
        T_MFPT = P_total_sJ/(1-P_off);
   end
end 


function [t,P] = solve_ode(glob_T,P_0,delta,TT)
     % glob_T is the given matrix A; initial condition is P_0;
    tspan = 0:delta:TT; % Time interval
    [t,P] = ode45(@(t,P)glob_T* P,tspan, P_0); % Solve the differential equation
end
