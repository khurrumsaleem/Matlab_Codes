% IMPLEMENTATION OF THE RESPONSE MATRIX METHOD FOR FIXED SOURCE, ISOTROPIC 
% SCATTERING, ONE-GROUP OF ENERGY DISCRETE-ORDINATES PROBLEMS IN SLAB
% GEOMETRY. (REF: Silva et. al., 2013)

% INPUT: (PROBLEM SETUP)
%        N   = QUADRATURE ORDER
%        ZON = ZONE PARAMETERS
%        DOM = DOMAIN PARAMETERS
%        BC  = BOUNDARY CONDITIONS
%        TOL = TOLERANCE (STOP CRITERIA)

% OUTPUT: SCALAR_FLUX = SCALAR FLUX AT THE EDGES
%         ANGULAR_FLUX = ANGULAR FLUX AT THE EDGES
%         ITER = NUMBER OF ITERATIONS
%         TIME = CALCULATION TIME

function [SCALAR_FLUX, ANGULAR_FLUX, ITER, TIME] = MET_1D_RM(N, ZON, DOM, BC, TOL)

  % VARIABLES
  NR = length(DOM(1, :)); ntc = sum(DOM(2, :)); QUAD = GAUSS_QUAD(N);
  [R, P] = response_matrix(QUAD, ZON, DOM);
  err = 1; ITER = 0;

  % INICIALIZATION
  ANGULAR_FLUX = zeros(ntc + 1, N); SCALAR_FLUX = zeros(ntc + 1, 1);
  for m = 1: N
    if (m <= N/2 && BC(1) ~= -1)
      ANGULAR_FLUX(1, m) = BC(1);
    elseif (BC(2) ~= -1)
      ANGULAR_FLUX(ntc + 1, m) = BC(2);
    end
  end
  
  % ITERATIVE PROCESS
  tic;
  fprintf("ITER\t\tMAX_DEV(SFLUX)\n");
  while err > TOL
    err = 0; ITER = ITER + 1;
    
    % 1. LEFT -> RIGHT SWEEP
    IB = 0;
    for r = 1: NR
      nc = DOM(2, r);
      RM = R(:, :, r); PM = P(:, r);
      for i = 1: nc
        IB = IB + 1; IF = IB + 1;
        IN = [ANGULAR_FLUX(IB, 1: N / 2) ANGULAR_FLUX(IF, N / 2 + 1: N)]';
        OUT = RM * IN + PM; 
        ANGULAR_FLUX(IF, 1: N / 2) = OUT(1: N / 2);
        ANGULAR_FLUX(IB, N / 2 + 1: N) = OUT(N / 2 + 1: N);
      end
    end
    
    % IN CASE OF BOUNDARY CONDITION ON THE RIGHT 
    if (BC(2) == -1)
      for m = 1: N/2
        ANGULAR_FLUX(ntc + 1, N/2 + m) = ANGULAR_FLUX(ntc + 1, m);
      end
    end
    
    % 2. RIGHT -> LEFT SWEEP
    IB = ntc + 1;
    for r = NR: -1: 1
      nc = DOM(2, r);
      RM = R(:, :, r); PM = P(:, r);
      for i = 1: nc
        IB = IB - 1; IF = IB + 1;
        IN = [ANGULAR_FLUX(IB, 1: N / 2) ANGULAR_FLUX(IF, N / 2 + 1: N)]';
        OUT = RM * IN + PM; 
        ANGULAR_FLUX(IF, 1: N / 2) = OUT(1: N / 2);
        ANGULAR_FLUX(IB, N / 2 + 1: N) = OUT(N / 2 + 1: N);
      end
    end
    
    % IN CASE OF BOUNDARY CONDITION ON THE LEFT
    if (BC(1) == -1)
      for m = 1: N/2
        ANGULAR_FLUX(1, m) = ANGULAR_FLUX(1, N/2 + m);
      end
    end
    
    % 3. STOP CRITERIA
    IB = 0;
    for r = 1: NR
      nc = DOM(2, r);
      for i = 1: nc
        IB = IB + 1; IF = IB + 1;
        faux = SCALAR_FLUX(IB);
        sflux = 0; sflux2 = 0;
        for m = 1: N
          w = QUAD(m,2); 
          sflux = sflux + ANGULAR_FLUX(IB, m) * w;
          sflux2 = sflux2 + ANGULAR_FLUX(IF, m) * w;
        end
        sflux = 0.5 * sflux; sflux2 = 0.5 * sflux2;
        SCALAR_FLUX(IB) = sflux; SCALAR_FLUX(IF) = sflux2;
        if (abs(1 - faux/sflux) > err)
          err = abs(1 - faux/sflux);
        end
      end
    end
    
    fprintf("%d\t\t\t%.5e\n", ITER, err);
    
  end
  
  TIME = toc;

end
