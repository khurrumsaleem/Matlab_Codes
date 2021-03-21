% PROBLEM NO. 1 (Larsen, 1986)

  % QUADRATURE ORDER
  N = 2;

  % NUMBER OF ZONES
  NZ = 1;  

  ZON = [1.00;    % SIGMA_T
         0.97];   % SIGMA_S

  % NUMBER OF REGIONS
  NR = 1;   

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