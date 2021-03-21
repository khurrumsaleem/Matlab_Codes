% PROBLEM TEST NO. 1 (Semi-infinite problem)

  % QUADRATURE ORDER
  N = 2;

  % NUMBER OF ZONES
  NZ = 1;  

  ZON = [1.0;    % SIGMA_T
         0.5];   % SIGMA_S

  % NUMBER OF REGIONS
  NR = 1;   

  DOM = [1;      % LENGHT
         100;    % NUMBER OF NODES
          1;     % ZONE INDEX
          0];    % INTERNAL SOURCE

  % BOUNDARY CONDITIONS
  BC = [ 1;        % LEFT
        -1];       % RIGHT (REFLECTIVE)
    
  % TOLERANCE
  TOL = 1e-06;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%