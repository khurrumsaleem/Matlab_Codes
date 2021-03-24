% PROBLEM TEST NO. 1 (Symmetric problem without scattering)

  % QUADRATURE ORDER
  N = 2;

  % ZONES
  ZON = [1.0;    % SIGMA_T
         0.0];   % SIGMA_S

  % DOMAIN
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
