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
      
      RYIN = zeros(M, M); RYOUT = zeros(M, M); SY = zeros(M, 1);
      FY0 = zeros(M, M); FY1IN = zeros(M, M); FY1OUT = zeros(M, M);
      AY0 = zeros(M, M); AY1_PLUS = zeros(M, M); AY1_MINUS = zeros(M, M);
      
      SXX = zeros(M,1); XNABLA = zeros(M,M); XF0 = zeros(M,M); XF1 = zeros(M,M);
      FXX0 = zeros(M, M); FXX1IN = zeros(M, M); FXX1OUT = zeros(M, M);
      AXX0 = zeros(M, M); AXX1_PLUS = zeros(M, M); AXX1_MINUS = zeros(M, M);
    
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
            
            RYOUT(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            RYIN(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > M / 4 && m <= M / 2) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            RYOUT(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            RYIN(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > M / 2 && m <= 3 * M / 4) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            RYOUT(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            RYIN(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
          end
          
        elseif (m > 3 * M / 4 && m <= M) % ROWS
            
          for k = 1: M % COLUMNS
            RXOUT(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            RXIN(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            RYOUT(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            RYIN(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
          end
          
        end
        
        % SOURCE MATRICES
        SX(m) = Q / (st * (1 - c0));
        SY(m) = Q / (st * (1 - c0));
        SXX = 6 * m_theta * Q / (st * hy * st * (1 - c0));
        
        % SECONDARY MATRICES
        for k = 1: M % COLUMNS
            
          k_miu = QUAD(k, 1); k_theta = QUAD(k, 2); k_w = QUAD(k, 3);
          
          AX0(m,k) = - 0.25 * c0 * k_w * k_theta / (st * hy * (1 - c0));
          AX1_PLUS(m,k) = 0.5 * c0 * k_w * k_theta * (k_miu / (st * hx) + m_miu / (st * hx) + 0.5) / (st * hy * (1 - c0));
          AX1_MINUS(m,k) = 0.5 * c0 * k_w * k_theta * (k_miu / (st * hx) + m_miu / (st * hx) - 0.5) / (st * hy * (1 - c0));
          
          AY0(m,k) = - 0.25 * c0 * k_w * k_miu / (st * hx * (1 - c0));
          AY1_PLUS(m,k) = 0.5 * c0 * k_w * k_miu * (k_theta / (st * hy) + m_theta / (st * hy) + 0.5) / (st * hx * (1 - c0));
          AY1_MINUS(m,k) = 0.5 * c0 * k_w * k_miu * (k_theta / (st * hy) + m_theta / (st * hy) - 0.5) / (st * hx * (1 - c0));
          
          AXX0(m,k) = - 3 * 0.5 * c0 * k_w * k_theta * (0.5 + m_theta / (st*hy)) / (st * hy * (1 - c0)) ...
                      - 3 * 0.5 * c0 * k_w * k_theta * k_theta / (st * hy * st * hy * (1 - c0));
          AXX1_PLUS(m,k) = 3 * 0.5 * c0 * m_miu * k_w * k_theta * (1 + 4 * m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                           + 3 * c0 * m_miu * k_w * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                           + 3 * c0 * k_w * k_miu * k_theta * (0.5 + m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                           + 6 * c0 * k_w * k_miu * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                           + 3 * 0.5 * c0 * k_w * k_theta * (0.5 + 2 * m_theta / (st * hy)) / (st * hy * (1 - c0)) ...
                           + 3 * 0.5 * c0 * k_w * k_theta * k_theta / (st * hy * st * hy * (1 - c0));
          AX1_MINUS(m,k) = 3 * 0.5 * c0 * m_miu * k_w * k_theta * (1 + 4 * m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                           + 3 * c0 * m_miu * k_w * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                           + 3 * c0 * k_w * k_miu * k_theta * (0.5 + m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                           + 6 * c0 * k_w * k_miu * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                           - 3 * 0.5 * c0 * k_w * k_theta * (0.5 + 2 * m_theta / (st * hy)) / (st * hy * (1 - c0)) ...
                           - 3 * 0.5 * c0 * k_w * k_theta * k_theta / (st * hy * st * hy * (1 - c0));
          
          
          if (k == m)
              AX0(m,k) = AX0(m,k) - m_theta / (st * hy);
              AX1_PLUS(m,k) = AX1_PLUS(m,k) + 2 * m_theta * (m_miu / (st * hx) + 0.5) / (st * hy);
              AX1_MINUS(m,k) = AX1_MINUS(m,k) + 2 * m_theta * (m_miu / (st * hx) - 0.5) / (st * hy);
              
              AY0(m,k) = AY0(m,k) - m_miu / (st * hx);
              AY1_PLUS(m,k) = AY1_PLUS(m,k) + 2 * m_miu * (m_theta / (st * hy) + 0.5) / (st * hx);
              AY1_MINUS(m,k) = AY1_MINUS(m,k) + 2 * m_miu * (m_theta / (st * hy) - 0.5) / (st * hx);
              
              AXX0(m,k) = AXX0(m,k) - 3 * m_theta * (1 + 2 * m_theta / (st * hy)) / (st * hy);
              AXX1_PLUS(m,k) = AXX1_PLUS(m,k) ...
                               + 6 * m_miu * m_theta * (1 + 4 * m_theta / (st * hy))/ (st * hx * st * hy) ...
                               + 3 * m_theta * (1 + 2 * m_theta / (st * hy)) / (st * hy);
              AXX1_MINUS(m,k) = AXX1_MINUS(m,k) ...
                               + 6 * m_miu * m_theta * (1 + 4 * m_theta / (st * hy))/ (st * hx * st * hy) ...
                               - 3 * m_theta * (1 + 2 * m_theta / (st * hy)) / (st * hy);
          end
          
          if (k <= M / 4)
            FX0(m,k) = AX0(m,k);
            FX1IN(m,k) = AX1_PLUS(m,k);
            FX1OUT(m,k) = AX1_MINUS(m,k);
            
            FY0(m,k) = AY0(m,k);
            FY1IN(m,k) = AY1_PLUS(m,k);
            FY1OUT(m,k) = AY1_MINUS(m,k);
            
            FXX0(m,k) = AXX0(m,k);
            FXX1IN(m,k) = AXX1_PLUS(m,k);
            FXX1OUT(m,k) = AXX1_MINUS(m,k);
          elseif (k > M / 4 && k <= M / 2)
            FX0(m,k) = AX0(m,k);
            FX1IN(m,k) = AX1_MINUS(m,k);
            FX1OUT(m,k) = AX1_PLUS(m,k);
            
            FY0(m,k) = - AY0(m,k);
            FY1IN(m,k) = - AY1_PLUS(m,k);
            FY1OUT(m,k) = - AY1_MINUS(m,k);
            
            FXX0(m,k) = AXX0(m,k);
            FXX1IN(m,k) = AXX1_MINUS(m,k);
            FXX1OUT(m,k) = AXX1_PLUS(m,k);
          elseif (k > M / 2 && k <= 3 * M / 4)
            FX0(m,k) = - AX0(m,k);
            FX1IN(m,k) = - AX1_MINUS(m,k);
            FX1OUT(m,k) = - AX1_PLUS(m,k);
            
            FY0(m,k) = - AY0(m,k);
            FY1IN(m,k) = - AY1_MINUS(m,k);
            FY1OUT(m,k) = - AY1_PLUS(m,k);
            
            FXX0(m,k) = - AXX0(m,k);
            FXX1IN(m,k) = - AXX1_MINUS(m,k);
            FXX1OUT(m,k) = - AXX1_PLUS(m,k);
          elseif (k > 3 * M / 4 && k <= M)
            FX0(m,k) = - AX0(m,k);
            FX1IN(m,k) = - AX1_PLUS(m,k);
            FX1OUT(m,k) = - AX1_MINUS(m,k);
            
            FY0(m,k) = AY0(m,k);
            FY1IN(m,k) = AY1_MINUS(m,k);
            FY1OUT(m,k) = AY1_PLUS(m,k);
            
            FXX0(m,k) = - AXX0(m,k);
            FXX1IN(m,k) = - AXX1_PLUS(m,k);
            FXX1OUT(m,k) = - AXX1_MINUS(m,k);
          end
          
        end
        
        % LAMBDA AND GAMA
        for k = 1: M
          XNABLA(m,k) = 6 * m_theta * xvects(m,k,z) / (hy * m_miu);
          if (xvals(m,z) == xvals(k,z))
            XF0(m,k) = 0;
            XF1(m,k) = 1;
          else
            XF0(m,k) = xvals(m,z) * xvals(k,z) / (st * (xvals(m,z) - xvals(k,z)));
            XF1(m,k) = 0;
          end
        end
        
      end
          
      ZERO = zeros(M, M);
      IDEN = eye(M, M);
      RXIN_INV = inv(RXIN);
      RYIN_INV = inv(RYIN);
      
      % LAMBDA AND GAMA
      LAMBDA_XIN = zeros(M, M); GAMA_XIN = zeros(M, M);
      LAMBDA_XOUT = zeros(M, M); GAMA_XOUT = zeros(M, M);
      XINV = inv(xvects(:,:,z));
      XB = XINV * XNABLA;
      xlanda = xvects(:,:,z)*(XB.*XF0);
      xgama = xvects(:,:,z)*(XB.*XF1);
      for m = 1: M
        for k = 1: M
          if (m <= M/4)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = - xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
          elseif (m > M/4 && m <= M/2)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = - xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
          elseif (m > M/2 && m <= 3*M/4)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = - xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
          elseif (m > 3*M/4 && m <= M)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = - xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
          end
        end
      end
      
      % EQUATION 1
      MX1 = RXOUT * RXIN_INV;
      MX2 = FX0 - (RXOUT * RXIN_INV) * FX0;
      MX4 = FX1OUT - (RXOUT * RXIN_INV) * FX1IN;
      MX5 = (IDEN - RXOUT * RXIN_INV) * SX;
      
      % EQUATION 2
      MY1 = FY0 - (RYOUT * RYIN_INV) * FY0;
      MY2 = RYOUT * RYIN_INV;
      MY3 = FY1OUT - (RYOUT * RYIN_INV) * FY1IN;
      MY5 = (IDEN - RYOUT * RYIN_INV) * SY;
      
      % EQUATION 3
      MXX1 = (- RXOUT * RXIN_INV * (LAMBDA_XIN + GAMA_XIN) + (LAMBDA_XOUT + GAMA_XOUT)) * RXIN_INV;
      MXX2 = - RXOUT * RXIN_INV * (FXX0 - (LAMBDA_XIN + GAMA_XIN) * RXIN_INV * FX0) ...
             + FXX0 - (LAMBDA_XOUT + GAMA_XOUT) * RXIN_INV * FX0;
      MXX3 = RXOUT * RXIN_INV;
      MXX4 = - RXOUT * RXIN_INV * (FXX1IN - (LAMBDA_XIN + GAMA_XIN) * RXIN_INV * FX1IN) ...
             + FXX1OUT - (LAMBDA_XOUT + GAMA_XOUT) * RXIN_INV * FX1IN;
      MXX5 = - RXOUT * RXIN_INV * (FXX0 - (LAMBDA_XIN + GAMA_XIN) * RXIN_INV * SX) ...
             + SXX - (LAMBDA_XOUT + GAMA_XOUT) * RXIN_INV * SX;
      
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