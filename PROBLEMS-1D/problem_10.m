% PROBLEM NO. 1 (Larsen, 1986)
 
  % ZONES
  ZON = [1.00;    % SIGMA_T
         0.97];   % SIGMA_S

  % DOMAIN
  DOM = [100;     % LENGHT
         100;     % NUMBER OF NODES
           1;     % ZONE INDEX
         0.0];    % INTERNAL SOURCE

  % BOUNDARY CONDITIONS
  BC = [1;        % LEFT
        0];       % RIGHT
    
  % TOLERANCE
  TOL = 1e-06;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
