% PROBLEM TEST NO. 1 (Symmetric problem without scattering)

  % QUADRATURE ORDER
  N = 2;

  % NUMBER OF ZONES
  NZ = 1;  

  ZON = [1.0;    % SIGMA_T
         0.0];   % SIGMA_S

  % NUMBER OF REGIONS
  NR = 1;   

  DOM = [1;      % LENGHT
         100;    % NUMBER OF NODES
          1;     % ZONE INDEX
          1];    % INTERNAL SOURCE

  % BOUNDARY CONDITIONS
  BC = [0;        % LEFT
        0];       % RIGHT
    
  % TOLERANCE
  TOL = 1e-06;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%