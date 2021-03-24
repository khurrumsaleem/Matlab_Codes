% PROBLEM TEST NO. 1 (Semi-infinite problem)

  % QUADRATURE ORDER
  N = 2;

  % ZONES
  ZON = [1.0;    % SIGMA_T
         0.5];   % SIGMA_S

  % DOMAIN
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
