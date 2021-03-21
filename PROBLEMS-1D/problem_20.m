% PROBLEM NO. 2 (Ortiz et. al., 2019)

  % QUADRATURE ORDER
  N = 8;

  % NUMBER OF ZONES
  NZ = 2;  

  ZON = [1.0, 0.6;    % SIGMA_T
         0.9, 0.4];   % SIGMA_S

  % NUMBER OF REGIONS
  NR = 3;   

  DOM = [20, 50, 30;     % LENGHT
         20, 50, 30;     % NUMBER OF NODES
          1,  2,  1;     % ZONE INDEX
          0,  0,  0];    % INTERNAL SOURCE

  % BOUNDARY CONDITIONS
  BC = [1;        % LEFT
        0];       % RIGHT
    
  % TOLERANCE
  TOL = 1e-06;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%