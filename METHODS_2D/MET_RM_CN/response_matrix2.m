% CONSTRUCTION OF THE RESPONSE MATRIX R AND THE SOURCE VECTOR P FOR THE
% ITERATIVE PROCESS

% OUTPUT: R = RESPONSE MATRIX FOR THE ITERATIVE PROCESS
%         P = SOURCE VECTOR FOR THE ITERATIVE PROCESS
%         RM = RESPONSE MATRIZ FOR THE RECONSTRUCTION OF THE AVERAGE
%              ANGULAR FLUXES IN EACH NODE
%         LM = LEAKAGE MATRIX FOR THE RECONTRUCTION OF THE AVERAGE ANGULAR
%              FLUXES IN EACH NODE
%         PM = SOURCE VECTOR FOR THE RECONSTRUCTION OF THE AVERAGE ANGULAR
%              FLUXES IN EACH NODE

function [R, P, F0, F1, SP] = response_matrix2(N, ZON, XDOM, YDOM, ZMAP, QMAP)
% SHORTHAND
% R(m, k, ry, rx)   P(m, ry, rx)
  
  % INITIALIZE VARIABLES
  NRX = length(XDOM(1, :)); NRY = length(YDOM(1, :));
  [QUAD, chi] = LQN_2D(N);
  [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
  
  % RESPONSE MATRIX AND SOURCE VECTOR
  M = N * (N + 2) / 2;
  R = zeros(2 * M, 2 * M, NRY, NRX); P = zeros(2 * M, NRY, NRX);
  F0 = zeros(M, M, NRY, NRX); F1 = zeros(M, M, NRY, NRX);
  SP = zeros(M, NRY, NRX);
  
  % AUXILIARY MATRICES
  XA = zeros(M, M); YA = zeros(M, M); XB = zeros(M, M); YB = zeros(M, M);
  XE = zeros(M, M); YE = zeros(M, M); I = eye(M, M);    HX = zeros(M, 1); HY = zeros(M, 1);
  
  
  % DEFINE EACH MATRIX PER REGION
  for ry = 1: NRY
    for rx = 1: NRX
    
      % FILL AUXILIARY MATRICES
      lenx = XDOM(1, rx); ntcx = XDOM(2, rx); hx = lenx / ntcx;
      leny = YDOM(1, ry); ntcy = YDOM(2, ry); hy = leny / ntcy;
      z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
      Q = QMAP(ry, rx);
      for m = 1: M
          
        % MATRICES FOR AVERAGE ANGULAR FLUX
        m_miu = QUAD(m,1); m_theta = QUAD(m,2);
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
        
        % XA, YA, XE, YE MATRICES
        if (m <= M / 4) % ROWS
          for k = 1: M  % COLUMNS
            XA(m, k) =  xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            XE(m, k) =  xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > M / 4 && m <= M / 2) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            XE(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > M / 2 && m <= 3 * M / 4) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            XE(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > 3 * M / 4 && m <= M) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            XE(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YE(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            
          end
        end
        
        % XB, YB MATRICES
        for k = 1: M % COLUMNS
          miu = QUAD(k, 1); theta = QUAD(k, 2); w = QUAD(k, 3);
          if (k <= M / 4)
            XB(m, k) = - 0.25 * c0 * theta * w / (1 - c0);
            if (k == m), XB(m, k) = XB(m,k) - theta; end
         
            YB(m, k) = - 0.25 * c0 * miu * w / (1 - c0);
            if (k == m), YB(m, k) = YB(m,k) - miu ; end
             
          elseif (k > M / 4 && k <= M / 2)
            XB(m, k) = - 0.25 * c0 * theta * w / ( (1 - c0));
            if (k == m), XB(m, k) = XB(m,k) - theta; end
             
            YB(m, k) = 0.25 * c0 * miu * w / ((1 - c0));
            if (k == m), YB(m, k) = YB(m,k) + miu; end
              
          elseif (k > M / 2 && k <= 3 * M / 4)
            XB(m, k) = 0.25 * c0 * theta * w / ((1 - c0));
            if (k == m), XB(m, k) = XB(m,k) + theta; end
              
            YB(m, k) = 0.25 * c0 * miu * w / ((1 - c0));
            if (k == m), YB(m, k) = YB(m,k) + miu; end
              
          elseif (k > 3 * M / 4 && k <= M)
            XB(m, k) = 0.25 * c0 * theta * w /(1 - c0);
            if (k == m), XB(m, k) = XB(m,k) + theta; end
              
            YB(m, k) = - 0.25 * c0 * miu * w / ((1 - c0));
            if (k == m), YB(m, k) = YB(m,k) - miu; end
              
          end
          
          
        end
        
        % H COLUMN VECTOR
        HX(m) = Q / (st * (1 - c0));
        HY(m) = Q / (st * (1 - c0));
        
        
      end
      
      
    
      XE_INV = inv(XE); 
      YE_INV = inv(YE);
      
      AUX = [           I.*(st*hx)                 , (- (XB - XA * (XE_INV * XB)));
             (- (YB - YA * (YE_INV * YB))),              I.*(st*hy)              ];
      
      AUX_INV = inv(AUX);
      
      AUX2 = [        (XA * XE_INV).*(st*hx)        , (- (XB - XA * (XE_INV * XB)));
              (- (YB - YA * (YE_INV * YB))),       (YA * YE_INV).*(st*hy)       ];
      
      % CALCULATE THE RESPONSE MATRIX
      R(:, :, ry, rx) = AUX_INV * AUX2;
      
      AUX3 = [(I - XA * XE_INV),    zeros(M, M)    ;
                zeros(M, M)    , (I - YA * YE_INV)];
      
      S = [HX; HY];
      
      % CALCULATE THE SOURCE VECTOR
      P(:, ry, rx) = AUX_INV * (AUX3 * S);
    
    end
  end
  
end