% CONSTRUCTION OF THE RESPONSE MATRIX R AND THE SOURCE VECTOR P FOR THE
% ITERATIVE PROCESS

% OUTPUT: R = RESPONSE MATRIX FOR THE ITERATIVE PROCESS
%         S = SOURCE VECTOR FOR THE ITERATIVE PROCESS
%         RM = RESPONSE MATRIZ FOR THE RECONSTRUCTION OF THE AVERAGE
%              ANGULAR FLUXES IN EACH NODE
%         LM = LEAKAGE MATRIX FOR THE RECONTRUCTION OF THE AVERAGE ANGULAR
%              FLUXES IN EACH NODE
%         PM = SOURCE VECTOR FOR THE RECONSTRUCTION OF THE AVERAGE ANGULAR
%              FLUXES IN EACH NODE

function [R, S] = response_matrix(N, ZON, XDOM, YDOM, ZMAP, QMAP)
% SHORTHAND
% 
  
  % INITIALIZE VARIABLES
  NRX = length(XDOM(1, :)); NRY = length(YDOM(1, :));
  [QUAD, chi] = LQN_2D(N);
  [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
  
  % RESPONSE MATRIX AND SOURCE VECTOR
  M = N * (N + 2) / 2;
  R = zeros(4 * M, 4 * M, NRY, NRX); S = zeros(4 * M, NRY, NRX);
  
  % DEFINE EACH MATRIX PER REGION
  for ry = 1: NRY
    for rx = 1: NRX
        
      % AUXILIARY MATRICES NEW
      RXIN = zeros(M, M); RXOUT = zeros(M, M); SX = zeros(M,1);
      FX0 = zeros(M, M); FX1IN = zeros(M, M); FX1OUT = zeros(M, M);
      AX0 = zeros(M, M); AX1_PLUS = zeros(M, M); AX1_MINUS = zeros(M, M);
    
      % AUXILIARY CONSTANTS
      lenx = XDOM(1, rx); ntcx = XDOM(2, rx); hx = lenx / ntcx;
      leny = YDOM(1, ry); ntcy = YDOM(2, ry); hy = leny / ntcy;
      z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
      Q = QMAP(ry, rx);
      
      % DIRECTION LOOP
      for m = 1: M
        
        m_miu = QUAD(m, 1); m_theta = QUAD(m, 2); m_w = QUAD(m, 3);
          
        % PRIMARY MATRICES
        if (m <= M / 4) % ROWS
            
          for k = 1: M  % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > M / 4 && m <= M / 2) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > M / 2 && m <= 3 * M / 4) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > 3 * M / 4 && m <= M) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
          end
          
        end
        
        % SOURCE MATRICES
        SX(m) = Q / (st * (1 - c0));
        
        % SECONDARY MATRICES
        for k = 1: M % COLUMNS
            
          k_miu = QUAD(k, 1); k_theta = QUAD(k, 2); k_w = QUAD(k, 3);
          
          AX0(m,k) = - 0.25 * c0 * k_w * k_theta / (st * hy * (1 - c0));
          AX1_PLUS(m,k) = 0.5 * c0 * k_w * k_theta * (k_miu / (st * hx) + m_miu / (st * hx) + 0.5) / (st * hy * (1 - c0));
          AX1_MINUS(m,k) = 0.5 * c0 * k_w * k_theta * (k_miu / (st * hx) + m_miu / (st * hx) - 0.5) / (st * hy * (1 - c0));
          if (k == m)
              AX0(m,k) = AX0(m,k) - m_theta / (st * hy);
              AX1_PLUS(m,k) = AX1_PLUS(m,k) + 2 * m_theta * (m_miu / (st * hx) + 0.5) / (st * hy);
              AX1_MINUS(m,k) = AX1_MINUS(m,k) + 2 * m_theta * (m_miu / (st * hx) - 0.5) / (st * hy);
          end
          
          if (k <= M / 4)
            FX0(m,k) = AX0(m,k);
            FX1IN(m,k) = AX1_PLUS(m,k);
            FX1OUT(m,k) = AX1_MINUS(m,k);
          elseif (k > M / 4 && k <= M / 2)
            FX0(m,k) = AX0(m,k);
            FX1IN(m,k) = AX1_MINUS(m,k);
            FX1OUT(m,k) = AX1_PLUS(m,k);
          elseif (k > M / 2 && k <= 3 * M / 4)
            FX0(m,k) = - AX0(m,k);
            FX1IN(m,k) = - AX1_MINUS(m,k);
            FX1OUT(m,k) = - AX1_PLUS(m,k);
          elseif (k > 3 * M / 4 && k <= M)
            FX0(m,k) = - AX0(m,k);
            FX1IN(m,k) = - AX1_PLUS(m,k);
            FX1OUT(m,k) = - AX1_MINUS(m,k);
          end
          
        end
        
      end
          
      ZERO = zeros(M, M);
      IDEN = eye(M, M);
      RXIN_INV = inv(RXIN);
      
      % EQUATION 1
      MX1 = RXOUT * RXIN_INV;
      MX2 = FX0 - (RXOUT * RXIN_INV) * FX0;
      MX4 = FX1OUT - (RXOUT * RXIN_INV) * FX1IN;
      MX5 = (IDEN - RXOUT * RXIN_INV) * SX;
      
      IN = [MX1,   -MX2,  ZERO,  -MX4;
            -MY1,  MY2,   -MY3,  ZERO;
            MXX1,  -MXX2, MXX3,  -MXX4;
            -MYY1, MYY2,  -MYY3, MYY4];
          
      OUT = [IDEN,  -MX2,  ZERO,  -MX4;
             -MY1,  IDEN,  -MY3,  ZERO;
             ZERO,  -MXX2, IDEN,  -MXX4;
             -MYY1, ZERO,  -MYY3, IDEN];
         
      OUT_INV = inv(OUT);
          
      % CALCULATE THE RESPONSE MATRIX
      R(:, :, ry, rx) = OUT_INV * IN;
      
      IND = [MX5; MY5; MXX5; MYY5];
      
      % CALCULATE THE SOURCE VECTOR
      S(:, ry, rx) = OUT_INV * IND;
    
    end
  end
  
end