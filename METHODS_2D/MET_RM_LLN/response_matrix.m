% CONSTRUCTION OF THE RESPONSE MATRIX R AND THE SOURCE VECTOR P FOR THE
% ITERATIVE PROCESS

function [R, S, F0, F1, SP] = response_matrix(N, ZON, XDOM, YDOM, ZMAP, QMAP)
  
  % INITIALIZE VARIABLES
  NRX = length(XDOM(1, :)); NRY = length(YDOM(1, :));
  [QUAD, chi] = LQN_2D(N);
  [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
  
  % RESPONSE MATRIX AND SOURCE VECTOR
  M = N * (N + 2) / 2;
  R = zeros(4 * M, 4 * M, NRY, NRX); S = zeros(4 * M, NRY, NRX);
  F0 = zeros(M, M, NRY, NRX); F1 = zeros(M, M, NRY, NRX);
  SP = zeros(M, NRY, NRX);
  
  % DEFINE EACH MATRIX PER REGION
  for ry = 1: NRY
    for rx = 1: NRX
        
      % AUXILIARY MATRICES
      RXIN = zeros(M, M); RXOUT = zeros(M, M); SX = zeros(M,1);
      FX0 = zeros(M, M); FX1IN = zeros(M, M); FX1OUT = zeros(M, M);
      AX0 = zeros(M, M); AX1_PLUS = zeros(M, M); AX1_MINUS = zeros(M, M);
      
      RYIN = zeros(M, M); RYOUT = zeros(M, M); SY = zeros(M, 1);
      FY0 = zeros(M, M); FY1IN = zeros(M, M); FY1OUT = zeros(M, M);
      AY0 = zeros(M, M); AY1_PLUS = zeros(M, M); AY1_MINUS = zeros(M, M);
      
      SXX = zeros(M,1); XNABLA = zeros(M,M); XF0 = zeros(M,M); XF1 = zeros(M,M);
      FXX0 = zeros(M, M); FXX1IN = zeros(M, M); FXX1OUT = zeros(M, M);
      AXX0 = zeros(M, M); AXX1_PLUS = zeros(M, M); AXX1_MINUS = zeros(M, M);
      
      SYY = zeros(M,1); YNABLA = zeros(M,M); YF0 = zeros(M,M); YF1 = zeros(M,M);
      FYY0 = zeros(M, M); FYY1IN = zeros(M, M); FYY1OUT = zeros(M, M);
      AYY0 = zeros(M, M); AYY1_PLUS = zeros(M, M); AYY1_MINUS = zeros(M, M);
    
      % AUXILIARY CONSTANTS
      lenx = XDOM(1, rx); ntcx = XDOM(2, rx); hx = lenx / ntcx;
      leny = YDOM(1, ry); ntcy = YDOM(2, ry); hy = leny / ntcy;
      z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
      Q = QMAP(ry, rx);
      
      % DIRECTION LOOP
      for m = 1: M
        
        m_miu = QUAD(m, 1); m_theta = QUAD(m, 2);
        
        % MATRICES FOR MEAN ANGULAR FLUX
        SP(m,ry,rx) = Q * (1 + c0 / (1 - c0)) / st;
        for k = 1: M
          k_miu = QUAD(k,1); k_theta = QUAD(k,2); k_w = QUAD(k,3);
          F0(m,k,ry,rx) = 0.25 * c0 * k_w * k_miu / (st * hx * (1-c0));
          F1(m,k,ry,rx) = 0.25 * c0 * k_w * k_theta / (st * hy * (1-c0));
          if (m == k)
            F0(m,k,ry,rx) = F0(m,k,ry,rx) + m_miu / (st * hx);
            F1(m,k,ry,rx) = F1(m,k,ry,rx) + m_theta / (st * hy);
          end
          if (k <= M/4)
            F0(m,k,ry,rx) = F0(m,k,ry,rx);
            F1(m,k,ry,rx) = F1(m,k,ry,rx);
          elseif (k > M/4 && k <= M/2)
            F0(m,k,ry,rx) = - F0(m,k,ry,rx);
            F1(m,k,ry,rx) = F1(m,k,ry,rx);
          elseif (k > M/2 && k <= 3*M/4)
            F0(m,k,ry,rx) = - F0(m,k,ry,rx);
            F1(m,k,ry,rx) = - F1(m,k,ry,rx);
          elseif (k > 3*M/4 && k <= M)
            F0(m,k,ry,rx) = F0(m,k,ry,rx);
            F1(m,k,ry,rx) = - F1(m,k,ry,rx);
          end
        end
          
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
        SXX(m) = 6 * m_theta * Q / (st * hy * st * (1 - c0));
        SYY(m) = 6 * m_miu * Q / (st * hx * st * (1 - c0));
        
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
          AXX1_MINUS(m,k) = 3 * 0.5 * c0 * m_miu * k_w * k_theta * (1 + 4 * m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                            + 3 * c0 * m_miu * k_w * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                            + 3 * c0 * k_w * k_miu * k_theta * (0.5 + m_theta / (st * hy)) / (st * hx * st * hy * (1 - c0)) ...
                            + 6 * c0 * k_w * k_miu * k_theta * k_theta / (st * hx * st * hy * st * hy * (1 - c0)) ...
                            - 3 * 0.5 * c0 * k_w * k_theta * (0.5 + 2 * m_theta / (st * hy)) / (st * hy * (1 - c0)) ...
                            - 3 * 0.5 * c0 * k_w * k_theta * k_theta / (st * hy * st * hy * (1 - c0));
                       
          AYY0(m,k) = - 3 * 0.5 * c0 * k_w * k_miu * (0.5 + m_miu / (st*hx)) / (st * hx * (1 - c0)) ...
                      - 3 * 0.5 * c0 * k_w * k_miu * k_miu / (st * hx * st * hx * (1 - c0));
          AYY1_PLUS(m,k) = 3 * 0.5 * c0 * m_theta * k_w * k_miu * (1 + 4 * m_miu / (st * hx)) / (st * hx * st * hy * (1 - c0)) ...
                           + 3 * c0 * m_theta * k_w * k_miu * k_miu / (st * hx * st * hx * st * hy * (1 - c0)) ...
                           + 3 * c0 * k_w * k_miu * k_theta * (0.5 + m_miu / (st * hx)) / (st * hx * st * hy * (1 - c0)) ...
                           + 6 * c0 * k_w * k_miu * k_miu * k_theta / (st * hx * st * hx * st * hy * (1 - c0)) ...
                           + 3 * 0.5 * c0 * k_w * k_miu * (0.5 + 2 * m_miu / (st * hx)) / (st * hx * (1 - c0)) ...
                           + 3 * 0.5 * c0 * k_w * k_miu * k_miu / (st * hx * st * hx * (1 - c0));
          AYY1_MINUS(m,k) = 3 * 0.5 * c0 * m_theta * k_w * k_miu * (1 + 4 * m_miu / (st * hx)) / (st * hx * st * hy * (1 - c0)) ...
                            + 3 * c0 * m_theta * k_w * k_miu * k_miu / (st * hx * st * hx * st * hy * (1 - c0)) ...
                            + 3 * c0 * k_w * k_miu * k_theta * (0.5 + m_miu / (st * hx)) / (st * hx * st * hy * (1 - c0)) ...
                            + 6 * c0 * k_w * k_miu * k_miu * k_theta / (st * hx * st * hx * st * hy * (1 - c0)) ...
                            - 3 * 0.5 * c0 * k_w * k_miu * (0.5 + 2 * m_miu / (st * hx)) / (st * hx * (1 - c0)) ...
                            - 3 * 0.5 * c0 * k_w * k_miu * k_miu / (st * hx * st * hx * (1 - c0));
          
          
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
                           
              AYY0(m,k) = AYY0(m,k) - 3 * m_miu * (1 + 2 * m_miu / (st * hx)) / (st * hx);
              AYY1_PLUS(m,k) = AYY1_PLUS(m,k) ...
                               + 6 * m_miu * m_theta * (1 + 4 * m_miu / (st * hx))/ (st * hx * st * hy) ...
                               + 3 * m_miu * (1 + 2 * m_miu / (st * hx)) / (st * hx);
              AYY1_MINUS(m,k) = AYY1_MINUS(m,k) ...
                               + 6 * m_miu * m_theta * (1 + 4 * m_miu / (st * hx))/ (st * hx * st * hy) ...
                               - 3 * m_miu * (1 + 2 * m_miu / (st * hx)) / (st * hx);
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
            
            FYY0(m,k) = AYY0(m,k);
            FYY1IN(m,k) = AYY1_PLUS(m,k);
            FYY1OUT(m,k) = AYY1_MINUS(m,k);
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
            
            FYY0(m,k) = - AYY0(m,k);
            FYY1IN(m,k) = - AYY1_PLUS(m,k);
            FYY1OUT(m,k) = - AYY1_MINUS(m,k);
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
            
            FYY0(m,k) = - AYY0(m,k);
            FYY1IN(m,k) = - AYY1_MINUS(m,k);
            FYY1OUT(m,k) = - AYY1_PLUS(m,k);
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
            
            FYY0(m,k) = AYY0(m,k);
            FYY1IN(m,k) = AYY1_MINUS(m,k);
            FYY1OUT(m,k) = AYY1_PLUS(m,k);
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
          
          YNABLA(m,k) = 6 * m_miu * yvects(m,k,z) / (hx * m_theta);
          if (yvals(m,z) == yvals(k,z))
            YF0(m,k) = 0;
            YF1(m,k) = 1;
          else
            YF0(m,k) = yvals(m,z) * yvals(k,z) / (st * (yvals(m,z) - yvals(k,z)));
            YF1(m,k) = 0;
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
      LAMBDA_YIN = zeros(M, M); GAMA_YIN = zeros(M, M);
      LAMBDA_YOUT = zeros(M, M); GAMA_YOUT = zeros(M, M);
      
      XINV = inv(xvects(:,:,z));
      XB = XINV * XNABLA;
      xlanda = xvects(:,:,z)*(XB.*XF0);
      xgama = xvects(:,:,z)*(XB.*XF1);
      
      YINV = inv(yvects(:,:,z));
      YB = YINV * YNABLA;
      ylanda = yvects(:,:,z)*(YB.*YF0);
      ygama = yvects(:,:,z)*(YB.*YF1);
      for m = 1: M
        for k = 1: M
          if (m <= M/4)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = - 0.5 * hx * xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = 0.5 * hx * xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            
            LAMBDA_YIN(m,k) = ylanda(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            LAMBDA_YOUT(m,k) = ylanda(m,k) * exp(0.5 * st * hy / yvals(k,z));
            GAMA_YIN(m,k) = - 0.5 * hy * ygama(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            GAMA_YOUT(m,k) = 0.5 * hy * ygama(m,k) * exp(0.5 * st * hy / yvals(k,z));
          elseif (m > M/4 && m <= M/2)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = 0.5 * hx * xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = - 0.5 * hx * xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            
            LAMBDA_YIN(m,k) = ylanda(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            LAMBDA_YOUT(m,k) = ylanda(m,k) * exp(0.5 * st * hy / yvals(k,z));
            GAMA_YIN(m,k) = - 0.5 * hy * ygama(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            GAMA_YOUT(m,k) = 0.5 * hy * ygama(m,k) * exp(0.5 * st * hy / yvals(k,z));
          elseif (m > M/2 && m <= 3*M/4)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = 0.5 * hx * xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = - 0.5 * hx * xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            
            LAMBDA_YIN(m,k) = ylanda(m,k) * exp(0.5 * st * hy / yvals(k,z));
            LAMBDA_YOUT(m,k) = ylanda(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            GAMA_YIN(m,k) = 0.5 * hy * ygama(m,k) * exp(0.5 * st * hy / yvals(k,z));
            GAMA_YOUT(m,k) = - 0.5 * hy * ygama(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
          elseif (m > 3*M/4 && m <= M)
            LAMBDA_XIN(m,k) = xlanda(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            LAMBDA_XOUT(m,k) = xlanda(m,k) * exp(0.5 * st * hx / xvals(k,z));
            GAMA_XIN(m,k) = - 0.5 * hx * xgama(m,k) * exp(- 0.5 * st * hx / xvals(k,z));
            GAMA_XOUT(m,k) = 0.5 * hx * xgama(m,k) * exp(0.5 * st * hx / xvals(k,z));
            
            LAMBDA_YIN(m,k) = ylanda(m,k) * exp(0.5 * st * hy / yvals(k,z));
            LAMBDA_YOUT(m,k) = ylanda(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
            GAMA_YIN(m,k) = 0.5 * hy * ygama(m,k) * exp(0.5 * st * hy / yvals(k,z));
            GAMA_YOUT(m,k) = - 0.5 * hy * ygama(m,k) * exp(- 0.5 * st * hy / yvals(k,z));
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
         
      % EQUATION 4
      MYY1 = - RYOUT * RYIN_INV * (FYY0 - (LAMBDA_YIN + GAMA_YIN) * RYIN_INV * FY0) ...
             + FYY0 - (LAMBDA_YOUT + GAMA_YOUT) * RYIN_INV * FY0;
      MYY2 = (- RYOUT * RYIN_INV * (LAMBDA_YIN + GAMA_YIN) + (LAMBDA_YOUT + GAMA_YOUT)) * RYIN_INV;
      MYY3 = - RYOUT * RYIN_INV * (FYY1IN - (LAMBDA_YIN + GAMA_YIN) * RYIN_INV * FY1IN) ...
             + FYY1OUT - (LAMBDA_YOUT + GAMA_YOUT) * RYIN_INV * FY1IN;
      MYY4 = RYOUT * RYIN_INV;
      MYY5 = - RYOUT * RYIN_INV * (SYY - (LAMBDA_YIN + GAMA_YIN) * RYIN_INV * SY) ...
             + SYY - (LAMNDA_YOUT + GAMA_YOUT) * RYIN_INV * SY;
      
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