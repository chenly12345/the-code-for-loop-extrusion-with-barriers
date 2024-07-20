function [monomer_position,position_left_all,position_right_all, t_absorbed,t_off] = Gillespie_two_position_3D(params,N,k_diff,d_off,right_loading_site,barrier_site,epsilon)
    % This function is to simulate the time it takes for the cohesin subunit to reach the boundary or dissociate for the 3D case.
    % physical barrier  
    
    % initial value
    t              = 0;
    t_off          = 0;
    t_absorbed     = 0;
    position_left  = right_loading_site-1;
    position_right = right_loading_site;
    reaction_num   = 1:6;     % six reactions of cohesin position-switching
    k_barrier      = k_diff * exp(-epsilon);
    
    position_left_all  = position_left;
    position_right_all = position_right;
      
    %% Initializing chromatin and cohesin position
    monomer_position    = zeros(params.monomer_num, 3);  
    cohesin_position    = [position_left,position_right];  
    
    Rouse_matrix          = initializeConnectivityMatrix(params); % initialize rouse chain, harmonic potential between consecutive monomers
    interaction_matrix    = updateInteractionMatrix(params,cohesin_position,Rouse_matrix); % Calculate interaction matrix based on cohesin position and Rouse chain matrix
    
    rand_array = params.b.*randn(params.monomer_num,params.dimension); %generate a random array to initialize chromatin positions
        for ps_idx = 2:params.monomer_num                        
            monomer_position(ps_idx,:) = monomer_position(ps_idx-1,:) + rand_array(ps_idx,:);
        end
    
    randn_array = randn(params.monomer_num,params.dimension,1000).*params.factor;
        for index = 1:1000
            monomer_position = monomer_position -interaction_matrix*(monomer_position*params.dt) + randn_array(:,:,index);
        end
    
     i = 1;
     randn_array  = randn(params.monomer_num,params.dimension,500000).*params.factor;
    
    %% main simulation loop
    while position_left ~= 0 || position_right ~= N  
        % % The two subunits are not adjacent and neither is at the boundary
         % the physical barrier is adjacent to the left of the left subunit
       if position_left ~= 0 && position_right ~= position_left + 1 && position_right ~= N  && position_left == barrier_site + 1 
            rates = [k_barrier,k_diff,k_diff,k_diff,d_off,d_off];
            % the physical barrier is adjacent to the right of the left subunit but not the left of right subunit
        elseif  position_left ~= 0 && position_right ~= position_left + 1 && position_right ~= position_left + 2 &&  position_right ~= N && barrier_site == position_left + 1 % 不相邻，且都不在边界 ,障碍在左亚基右端 2
            rates = [k_diff,k_barrier,k_diff,k_diff,d_off,d_off];
            % the physical barrier is adjacent  to the right of the left subunit and the left of the right subunit
        elseif  position_left ~= 0 && position_right ~= position_left + 1 && position_right == position_left + 2 &&  position_right ~= N && barrier_site == position_left + 1 % 不相邻，且都不在边界 ,障碍在左亚基右端 2
            rates = [k_diff,0,0,k_diff,d_off,d_off];
            % the physical barrier is adjacent  to the left of the right subunit
        elseif position_left ~= 0 && position_right ~= position_left + 1 && position_right ~= N  && position_right == barrier_site + 1 
            rates = [k_diff,k_diff,k_barrier,k_diff,d_off,d_off];
            % the physical barrier is adjacent to the right of the right subunit
        elseif position_left ~= 0 && position_right ~= position_left + 1 && position_right ~= N  && barrier_site == position_right + 1
            rates = [k_diff,k_diff,k_diff,k_barrier,d_off,d_off];
            % the physical barrier and subunits are not adjacent
        elseif position_left ~= 0 && position_right ~= position_left + 1 && position_right ~= N  && barrier_site ~= position_right + 1 && position_right ~= barrier_site + 1 && barrier_site ~= position_left + 1 && position_left ~= barrier_site + 1
            rates = [k_diff,k_diff,k_diff,k_diff,d_off,d_off];
    
            % % The two subunits are adjacent and neither is at the boundary
            % the physical barrier is adjacent to the left of the left subunit
        elseif position_left ~= 0 && position_right == position_left + 1 &&   position_right ~= N && position_left == barrier_site + 1 
            rates = [k_barrier,0,0,k_diff,d_off,d_off];
            % the physical barrier is adjacent  to the right of the right subunit
        elseif position_left ~= 0 && position_right == position_left + 1 &&   position_right ~= N && barrier_site == position_right + 1 
            rates = [k_diff,0,0,k_barrier,d_off,d_off];  
            % the physical barrier and subunits are not adjacent
        elseif position_left ~= 0 && position_right == position_left + 1 &&   position_right ~= N && position_left ~= barrier_site + 1 && barrier_site ~= position_right + 1
            rates = [k_diff,0,0,k_diff,d_off,d_off];
          
            % % The left subunit is on the left boundary, and the right subunit is neither adjacent nor on the right boundary
            % the physical barrier is adjacent to the left of the right subunit
        elseif position_left == 0 && position_right ~= 1 && position_right ~= N  && position_right == barrier_site + 1 
            rates = [0,0,k_barrier,k_diff,0,d_off];    
            % the physical barrier is adjacent  to the right of the right subunit
        elseif position_left == 0 && position_right ~= 1 && position_right ~= N && barrier_site == position_right + 1
            rates = [0,0,k_diff,k_barrier,0,d_off];    
            % the physical barrier and subunits are not adjacent
        elseif position_left == 0 && position_right ~= 1 && position_right ~= N  && position_right ~= barrier_site + 1  && barrier_site ~= position_right + 1
            rates = [0,0,k_diff,k_diff,0,d_off];
           
            %  % The left subunit is on the left boundary, and the right subunit is adjacent
            % the physical barrier is adjacent to the right of the right subunit
        elseif position_left == 0 && position_right == 1 && barrier_site == position_right + 1 
            rates = [0,0,0,k_barrier,0,d_off];
            % the physical barrier and subunits are not adjacent
        elseif position_left == 0 && position_right == 1 && barrier_site ~= position_right + 1 
            rates = [0,0,0,k_diff,0,d_off];
         
            % % The left subunit is not at the boundary, not adjacent to the right subunit, and the right subunit is at the boundary
            % the physical barrier is adjacent  to the left of the left subunit
        elseif position_left ~= 0 && position_left ~= N-1 && position_right == N && position_left == barrier_site + 1 
            rates = [k_barrier,k_diff,0,0,d_off,0];
            % the physical barrier is adjacent  to the right of the left subunit
        elseif position_left ~= 0 && position_left ~= N-1 && position_right == N && barrier_site == position_left + 1    
            rates = [k_diff,k_barrier,0,0,d_off,0];
            % the physical barrier and subunits are not adjacent
        elseif position_left ~= 0 && position_left ~= N-1 && position_right == N && position_left ~= barrier_site + 1  && barrier_site ~= position_left + 1 
            rates = [k_diff,k_diff,0,0,d_off,0];
    
            % % The left subunit is not at the boundary, adjacent to the right subunit, and the right subunit is at the boundary
            % the physical barrier is adjacent  to the left of the left subunit
        elseif  position_left == N-1 && position_right == N  && position_left == barrier_site + 1 
           rates = [k_barrier,0,0,k_diff,d_off,0];
            % the physical barrier and subunits are not adjacent
        elseif  position_left == N-1 && position_right == N  && position_left ~= barrier_site + 1 
           rates = [k_diff,0,0,0,d_off,0];
        end
        
        % Compute tau and mu using random variates
        a0  = sum(rates);
        r1  = rand;
        tau = (1/a0)*log(1/r1);
        r2  = rand;
    
        % Calculate reaction propensities
        for iter = 1:length(reaction_num)
            if (sum(rates(1:iter)) >= r2*a0)
                next_mu = reaction_num(iter);
                break;
            end
        end
        
            last_time    = t(end);
            current_time = t(end) + tau;
            steps        = (current_time - last_time)/params.dt;
    
    % %  update the total interaction matrix   
     
      for index = 1:steps
            % update the chromatin position
            monomer_position = monomer_position -interaction_matrix*(monomer_position*params.dt) + randn_array(:,:,i);
                i = i + 1; 
                if i == 500000 
                  randn_array  = randn(params.monomer_num,params.dimension,500000).*params.factor;
                     i = 1; 
                end
      end
    
        % update the position and time 
        if  next_mu == 1 && position_left == barrier_site + 1 
                position_left =  position_left -2;
            elseif next_mu == 1 && position_left ~= barrier_site + 1
                position_left =  position_left -1;
            elseif next_mu == 2  && barrier_site == position_left + 1   
                 position_left =  position_left +2;
            elseif next_mu == 2  && barrier_site ~= position_left + 1     
                position_left =  position_left +1;
            elseif next_mu == 3   && position_right == barrier_site + 1 
                 position_right = position_right -2;
            elseif next_mu == 3  && position_right ~= barrier_site + 1  
                 position_right = position_right -1;   
            elseif next_mu == 4  && barrier_site == position_right + 1   
                position_right =position_right +2;
            elseif next_mu == 4  && barrier_site ~= position_right + 1   
                position_right =position_right +1;
            else
                t_off =t(end)+tau;
                break ;
         end
    
        cohesin_position      = [position_left,position_right];
        interaction_matrix    = updateInteractionMatrix(params,cohesin_position,Rouse_matrix);
    
       % % update time 
            t                  = [t; t(end) + tau];
            position_left_all  = [position_left_all;position_left];
            position_right_all = [position_right_all;position_right];
    
        if position_left == 0 && position_right == N
            t_absorbed=t(end);
            break;
        end   
     end   %end the simulation
end