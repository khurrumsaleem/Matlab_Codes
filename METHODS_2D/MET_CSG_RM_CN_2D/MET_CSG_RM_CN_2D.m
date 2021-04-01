% IMPLEMENTATION OF THE COMPOSITE SPATIAL GRID - RESPONSE MATRIX - CONSTANT
% NODAL (RM-CN) METHOD FOR FIXED SOURCE, ISOTROPIC SCATTERING, ONE-GROUP OF
% ENERGY DISCRETE - ORDINATES PROBLEMS IN X,Y - GEOMETRY.

function [SCALAR_FLUX, ...  % SCALAR FLUX IN EACH COLUMN FOR X
          ANGULAR_FLUX, ... % ANGULAR FLUX IN EACH COLUMN FOR X
          X_ANG_FLUX, ...   % ANGULAR FLUX AT X EDGES
          Y_ANG_FLUX, ...   % ANGULAR FLUX AT Y EDGES
          ITER, ...         % NUMBER OF ITERATIONS
          TIME ...          % CPU TIME
          ] = MET_CSG_RM_CN_2D(N, ZON, XDOM, YDOM, ZMAP, QMAP, BC, TOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % VARIABLES
  M = N * (N + 2) / 2; I = sum(XDOM(2, :)); J = sum(YDOM(2, :));
  RX = length(XDOM(1, :)); RY = length(YDOM(1, :));
  [XR, XL, XP, YR, YL, YP, RM, LM, PM] = response_matrix(N, ZON, XDOM, YDOM, ZMAP, QMAP);
  [QUAD, ~] = LQN_2D(N);
  
  % INITIALIZATION
  SCALAR_FLUX = zeros(RY, I); ANGULAR_FLUX = zeros(M, RY, I);
  X_ANG_FLUX = zeros(M, RY, I + 1); Y_ANG_FLUX = zeros(M, J + 1, RX);
  for ry = 1: RY
    for m = 1: M
      if ( m <= M / 4 && BC(1) ~= -1)
        X_ANG_FLUX(m, ry, 1) = BC(1);
      elseif (m > M / 4 && m <= M / 2 && BC(3) ~= -1)
        X_ANG_FLUX(m, ry, I + 1) = BC(3);
      elseif (m > M / 2 && m <= 3 * M / 4 && BC(3) ~= -1)
        X_ANG_FLUX(m, ry, I + 1) = BC(3);
      elseif (m > 3 * M / 4 && m <= M && BC(1) ~= -1)
        X_ANG_FLUX(m, ry, 1) = BC(1);
      end
    end
  end
  for rx = 1: RX
    for m = 1: M
      if ( m <= M / 4 && BC(2) ~= -1)
        Y_ANG_FLUX(m, 1, rx) = BC(2);
      elseif (m > M / 4 && m <= M / 2 && BC(2) ~= -1)
        Y_ANG_FLUX(m, 1, rx) = BC(2);
      elseif (m > M / 2 && m <= 3 * M / 4 && BC(4) ~= -1)
        Y_ANG_FLUX(m, J + 1, rx) = BC(4);
      elseif (m > 3 * M / 4 && m <= M && BC(4) ~= -1)
        Y_ANG_FLUX(m, J + 1, rx) = BC(4);
      end
    end
  end
  ERR = 1; ITER = -1;
  
  % ITERATIVE PROCESS
  tic;
  fprintf("ITER\t\tMAX_DEV(SCALAR_FLUX)\n");
  while ERR > TOL
    ERR = 0; ITER = ITER + 1;
    X_FLUX = X_ANG_FLUX; Y_FLUX = Y_ANG_FLUX;
    
    
    count = 5;
    while count > 0
      count = count - 1;
      
      % IN CASE OF SYMMETRIC BOUNDARY CONDITIONS IN X
      if (BC(1) == -1 || BC(3) == -1)
        for ry = 1: RY
          for m = 1: M / 4
            if BC(1) == -1.0
              X_ANG_FLUX(m, ry, 1) = X_ANG_FLUX(M / 4 + m, ry, 1);
              X_ANG_FLUX(3 * M / 4 + m, ry, 1) = X_ANG_FLUX(M / 2 + m, ry, 1);
            end
            if BC(3) == -1.0
              X_ANG_FLUX(M / 4 + m, ry, I + 1) =  X_ANG_FLUX(m, ry, I + 1);
              X_ANG_FLUX(M / 2 + m, ry, I + 1) =  X_ANG_FLUX(3 * M / 4 + m, ry, I + 1);
            end
          end
        end
      end
      
      % 1. LEFT -> RIGHT SWEEP
      JB = 1;
      for ry = 1: RY
        ncy = YDOM(2, ry);
        JF = JB + ncy;
        IB = 0;
        for rx = 1: RX
          ncx = XDOM(2, rx);
          for i = 1: ncx
            IB = IB + 1; IF = IB + 1;
            
            XIN = [X_ANG_FLUX(1: M / 4, ry, IB);
                   X_ANG_FLUX(M / 4 + 1: M / 2, ry, IF);
                   X_ANG_FLUX(M / 2 + 1: 3 * M / 4, ry, IF);
                   X_ANG_FLUX(3 * M / 4 + 1: M, ry, IB)];
               
            YIN = [Y_FLUX(1: M / 4, JB, rx);
                   Y_FLUX(M / 4 + 1: M / 2, JB, rx);
                   Y_FLUX(M / 2 + 1: 3 * M / 4, JF, rx);
                   Y_FLUX(3 * M / 4 + 1: M, JF, rx)];
               
            YOUT = [Y_FLUX(1: M / 4, JF, rx); ...
                    Y_FLUX(M / 4 + 1: M / 2, JF, rx); ...
                    Y_FLUX(M / 2 + 1: 3 * M / 4, JB, rx); ...
                    Y_FLUX(3 * M / 4 + 1: M, JB, rx)];

            XOUT = XR(:, :, ry, rx) * XIN - XL(:, :, ry, rx) * (YOUT - YIN) ...
                   + XP(:, ry, rx);
            
            X_ANG_FLUX(1: M / 4, ry, IF) = XOUT(1: M / 4);
            X_ANG_FLUX(M / 4 + 1: M / 2, ry, IB) = XOUT(M / 4 + 1: M / 2);
            X_ANG_FLUX(M / 2 + 1: 3 * M / 4, ry, IB) = XOUT(M / 2 + 1: 3 * M / 4);
            X_ANG_FLUX(3 * M / 4 + 1: M, ry, IF) = XOUT(3 * M / 4 + 1: M);
          end
        end
        JB = JF;
      end
      
      % 2. LEFT <- RIGHT SWEEP
      JB = 1;
      for ry = 1: RY
        ncy = YDOM(2, ry);
        JF = JB + ncy;
        IB = I + 1;
        for rx = RX: -1: 1
          ncx = XDOM(2, rx);
          for i = 1: ncx
            IB = IB - 1; IF = IB + 1;
            
            XIN = [X_ANG_FLUX(1: M / 4, ry, IB);
                   X_ANG_FLUX(M / 4 + 1: M / 2, ry, IF);
                   X_ANG_FLUX(M / 2 + 1: 3 * M / 4, ry, IF);
                   X_ANG_FLUX(3 * M / 4 + 1: M, ry, IB)];
               
            YIN = [Y_FLUX(1: M / 4, JB, rx);
                   Y_FLUX(M / 4 + 1: M / 2, JB, rx);
                   Y_FLUX(M / 2 + 1: 3 * M / 4, JF, rx);
                   Y_FLUX(3 * M / 4 + 1: M, JF, rx)];
               
            YOUT = [Y_FLUX(1: M / 4, JF, rx); ...
                    Y_FLUX(M / 4 + 1: M / 2, JF, rx); ...
                    Y_FLUX(M / 2 + 1: 3 * M / 4, JB, rx); ...
                    Y_FLUX(3 * M / 4 + 1: M, JB, rx)];

            XOUT = XR(:, :, ry, rx) * XIN - XL(:, :, ry, rx) * (YOUT - YIN) ...
                   + XP(:, ry, rx);
            
            X_ANG_FLUX(1: M / 4, ry, IF) = XOUT(1: M / 4);
            X_ANG_FLUX(M / 4 + 1: M / 2, ry, IB) = XOUT(M / 4 + 1: M / 2);
            X_ANG_FLUX(M / 2 + 1: 3 * M / 4, ry, IB) = XOUT(M / 2 + 1: 3 * M / 4);
            X_ANG_FLUX(3 * M / 4 + 1: M, ry, IF) = XOUT(3 * M / 4 + 1: M);
          end
        end
        JB = JF;
      end
      
    end % END WHILE IN X
    
    
    count = 5;
    while count > 0
      count = count - 1;
      
      % IN CASE OF SYMMETRIC BOUNDARY CONDITIONS IN Y
      if(BC(2) == -1 || BC(4) == -1)
        for rx = 1: RX
          for m = 1: M / 4
            if BC(2) == -1.0
              Y_ANG_FLUX(m, 1, rx) = Y_ANG_FLUX(3 * M / 4 + m, 1, rx);
              Y_ANG_FLUX(M / 4 + m, 1, rx) = Y_ANG_FLUX(M / 2 + m, 1, rx);
            end
            if BC(4) == -1.0
              Y_ANG_FLUX(3 * M / 4 + m, J + 1, rx) = Y_ANG_FLUX(m, J + 1, rx);
              Y_ANG_FLUX(M / 2 + m, J + 1, rx) = Y_ANG_FLUX(M / 4 + m, J + 1, rx);
            end
          end
        end
      end
      
      % 3. BOTTOM -> TOP SWEEP
      IB = 1;
      for rx = 1: RX
        ncx = XDOM(2,rx);
        IF = IB + ncx;
        JB = 0;
        for ry = 1: RY
          ncy = YDOM(2, ry);
          for j = 1: ncy
            JB = JB + 1; JF = JB + 1;
            
            YIN = [Y_ANG_FLUX(1: M / 4, JB, rx);
                   Y_ANG_FLUX(M / 4 + 1: M / 2, JB, rx);
                   Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JF, rx);
                   Y_ANG_FLUX(3 * M / 4 + 1: M, JF, rx)];
          
            XIN = [X_FLUX(1: M / 4, ry, IB);
                   X_FLUX(M / 4 + 1: M / 2, ry, IF);
                   X_FLUX(M / 2 + 1: 3 * M / 4, ry, IF);
                   X_FLUX(3 * M / 4 + 1: M, ry, IB)];
             
            XOUT = [X_FLUX(1: M / 4, ry, IF);
                    X_FLUX(M / 4 + 1: M / 2, ry, IB);
                    X_FLUX(M / 2 + 1: 3 * M / 4, ry, IB);
                    X_FLUX(3 * M / 4 + 1: M, ry, IF)];
          
            YOUT = YR(:, :, ry, rx) * YIN - YL(:, :, ry, rx) * (XOUT - XIN) ...
                   + YP(:, ry, rx);
          
            Y_ANG_FLUX(1: M / 4, JF, rx) = YOUT(1: M / 4);
            Y_ANG_FLUX(M / 4 + 1: M / 2, JF, rx) = YOUT(M / 4 + 1: M / 2);
            Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JB, rx) = YOUT(M / 2 + 1: 3 * M / 4);
            Y_ANG_FLUX(3 * M / 4 + 1: M, JB, rx) = YOUT(3 * M / 4 + 1: M);
            
          end
        end
        IB = IF;
      end
      
      % 4. BOTTOM <- TOP SWEEP
      IB = 1;
      for rx = 1: RY
        ncx = XDOM(2,rx);
        IF = IB + ncx;
        JB = J + 1;
        for ry = RY: -1: 1
          ncy = YDOM(2, ry);
          for j = 1: ncy
            JB = JB - 1; JF = JB + 1;
            
            YIN = [Y_ANG_FLUX(1: M / 4, JB, rx);
                   Y_ANG_FLUX(M / 4 + 1: M / 2, JB, rx);
                   Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JF, rx);
                   Y_ANG_FLUX(3 * M / 4 + 1: M, JF, rx)];
          
            XIN = [X_FLUX(1: M / 4, ry, IB);
                   X_FLUX(M / 4 + 1: M / 2, ry, IF);
                   X_FLUX(M / 2 + 1: 3 * M / 4, ry, IF);
                   X_FLUX(3 * M / 4 + 1: M, ry, IB)];
             
            XOUT = [X_FLUX(1: M / 4, ry, IF);
                    X_FLUX(M / 4 + 1: M / 2, ry, IB);
                    X_FLUX(M / 2 + 1: 3 * M / 4, ry, IB);
                    X_FLUX(3 * M / 4 + 1: M, ry, IF)];
          
            YOUT = YR(:, :, ry, rx) * YIN - YL(:, :, ry, rx) * (XOUT - XIN) ...
                   + YP(:, ry, rx);
          
            Y_ANG_FLUX(1: M / 4, JF, rx) = YOUT(1: M / 4);
            Y_ANG_FLUX(M / 4 + 1: M / 2, JF, rx) = YOUT(M / 4 + 1: M / 2);
            Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JB, rx) = YOUT(M / 2 + 1: 3 * M / 4);
            Y_ANG_FLUX(3 * M / 4 + 1: M, JB, rx) = YOUT(3 * M / 4 + 1: M);
            
          end
        end
        IB = IF;
      end
      
    end % END WHILE IN Y
    
    
    % SCALAR AND ANGULAR FLUXES CALCULATION AND STOP CRITERIA
    JB = 1;
    for ry = 1: RY
      ncy = YDOM(2, ry);
      JF = JB + ncy;
      IB = 0;
      for rx = 1: RX
        ncx = XDOM(2, rx);
        for i = 1: ncx
          IB = IB + 1; IF = IB + 1;
            
          XIN = [X_ANG_FLUX(1: M / 4, ry, IB);
                 X_ANG_FLUX(M / 4 + 1: M / 2, ry, IF);
                 X_ANG_FLUX(M / 2 + 1: 3 * M / 4, ry, IF);
                 X_ANG_FLUX(3 * M / 4 + 1: M, ry, IB)];
               
          YIN = [Y_ANG_FLUX(1: M / 4, JB, rx);
                 Y_ANG_FLUX(M / 4 + 1: M / 2, JB, rx);
                 Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JF, rx);
                 Y_ANG_FLUX(3 * M / 4 + 1: M, JF, rx)];
               
          YOUT = [Y_ANG_FLUX(1: M / 4, JF, rx); ...
                  Y_ANG_FLUX(M / 4 + 1: M / 2, JF, rx); ...
                  Y_ANG_FLUX(M / 2 + 1: 3 * M / 4, JB, rx); ...
                  Y_ANG_FLUX(3 * M / 4 + 1: M, JB, rx)];

          OUTPUT = RM(:, :, ry, rx) * XIN - LM(:, :, ry, rx) * (YOUT - YIN) ...
                   + PM(:, ry, rx);
            
          ANGULAR_FLUX(1: M / 4, ry, IB) = OUTPUT(1: M / 4);
          ANGULAR_FLUX(M / 4 + 1: M / 2, ry, IB) = OUTPUT(M / 4 + 1: M / 2);
          ANGULAR_FLUX(M / 2 + 1: 3 * M / 4, ry, IB) = OUTPUT(M / 2 + 1: 3 * M / 4);
          ANGULAR_FLUX(3 * M / 4 + 1: M, ry, IB) = OUTPUT(3 * M / 4 + 1: M);
            
          flux = 0.25 * (OUTPUT' * QUAD(:, 3));
          flux0 = SCALAR_FLUX(ry, IB);
          SCALAR_FLUX(ry, IB) = flux;
          if abs(1 - flux0 / flux) > ERR, ERR = abs(1 - flux0 / flux); end
        end
      end
      JB = JF;
    end
    
    fprintf("%d\t\t\t%.5e\n", ITER, ERR);
    
  end
  
  TIME = toc;

end
