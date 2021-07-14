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

function [R, P, RM, LM, PM] = response_matrix(N, ZON, XDOM, YDOM, ZMAP, QMAP)
% SHORTHAND
% R(m, k, ry, rx)   P(m, ry, rx)
  
  % INITIALIZE VARIABLES
  NRX = length(XDOM(1, :)); NRY = length(YDOM(1, :));
  [QUAD, chi] = LQN_2D(N);
  [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
  
  % RESPONSE MATRIX AND SOURCE VECTOR
  M = N * (N + 2) / 2;
  R = zeros(2 * M, 2 * M, NRY, NRX); P = zeros(2 * M, NRY, NRX);
  RM = zeros(M, M, NRY, NRX); LM = zeros(M, M, NRY, NRX);
  PM = zeros(M, NRY, NRX);
  
  % AUXILIARY MATRICES
  XA = zeros(M, M); YA = zeros(M, M); XB = zeros(M, M); YB = zeros(M, M);
  XE = zeros(M, M); YE = zeros(M, M); I = eye(M, M);    H = zeros(M, 1);
  XM = zeros(M, M);
  
  % DEFINE EACH MATRIX PER REGION
  for ry = 1: NRY
    for rx = 1: NRX
    
      % FILL AUXILIARY MATRICES
      lenx = XDOM(1, rx); ntcx = XDOM(2, rx); hx = lenx / ntcx;
      leny = YDOM(1, ry); ntcy = YDOM(2, ry); hy = leny / ntcy;
      z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
      Q = QMAP(ry, rx);
      for m = 1: M
        
        % XA, YA, XE, YE MATRICES
        if (m <= M / 4) % ROWS
          for k = 1: M  % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            XE(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
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
            XB(m, k) = - 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m), XB(m, k) = XB(m,k) - theta / (st * hy); end
         
            YB(m, k) = - 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m), YB(m, k) = YB(m,k) - miu / (st * hx); end
             
          elseif (k > M / 4 && k <= M / 2)
            XB(m, k) = - 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m), XB(m, k) = XB(m,k) - theta / (st * hy); end
             
            YB(m, k) = 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m), YB(m, k) = YB(m,k) + miu / (st * hx); end
              
          elseif (k > M / 2 && k <= 3 * M / 4)
            XB(m, k) = 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m), XB(m, k) = XB(m,k) + theta / (st * hy); end
              
            YB(m, k) = 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m), YB(m, k) = YB(m,k) + miu / (st * hx); end
              
          elseif (k > 3 * M / 4 && k <= M)
            XB(m, k) = 0.25 * c0 * theta * w / (st * hy * (1 - c0));
            if (k == m), XB(m, k) = XB(m,k) + theta / (st * hy); end
              
            YB(m, k) = - 0.25 * c0 * miu * w / (st * hx * (1 - c0));
            if (k == m), YB(m, k) = YB(m,k) - miu / (st * hx); end
              
          end
          
          XM(m, k) = xvects(m, k, z) * xvals(k, z) * 2 * sinh(0.5 * st * hx / xvals(k, z)) / (st * hx);
          
        end
        
        % H COLUMN VECTOR
        H(m) = Q / (st * (1 - c0));
        
      end
    
      XE_INV = inv(XE); 
      YE_INV = inv(YE);
      
      % CALCULATE THE RESPONSE MATRIX FOR RECONSTRUCTION
      RM(:, :, ry, rx) = XM * XE_INV;
      
      % CALCULATE THE LEAKAGE MATRIX FOR RECONSTRUCTION
      LM(:, :, ry, rx) = XB - XM * (XE_INV * XB);
      
      % CALCULATE THE SOURCE VECTOR FRO RECONSTRUCTION
      PM(:, ry, rx) = (I - XM * XE_INV) * H;
      
      AUX = [           I                 , (- (XB - XA * (XE_INV * XB)));
             (- (YB - YA * (YE_INV * YB))),              I              ];
      
      AUX_INV = inv(AUX);
      
      AUX2 = [        (XA * XE_INV)        , (- (XB - XA * (XE_INV * XB)));
              (- (YB - YA * (YE_INV * YB))),       (YA * YE_INV)       ];
      
      % CALCULATE THE RESPONSE MATRIX
      R(:, :, ry, rx) = AUX_INV * AUX2;
      
      AUX3 = [(I - XA * XE_INV),    zeros(M, M)    ;
                zeros(M, M)    , (I - YA * YE_INV)];
      
      S = [H; H];
      
      % CALCULATE THE SOURCE VECTOR
      P(:, ry, rx) = AUX_INV * (AUX3 * S);
    
    end
  end
  
end