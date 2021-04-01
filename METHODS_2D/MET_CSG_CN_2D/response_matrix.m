% CONSTRUCTION OF THE RESPONSE MATRIX R, THE LEAKAGE MATRIX L, AND 
% THE SOURCE VECTOR P FOR THE ITERATIVE PROCESS

function [XR,... RESPONSE MATRIX IN X
          XL,... LEAKAGE MATRIX IN X
          XP,... SOURCE VECTOR IN X
          YR,... RESPONSE MATRIX IN Y
          YL,... LEAKAGE MATRIX IN Y
          YP,... SOURCE VECTOR IN Y
          RM,... RESPONSE MATRIX FOR RECONSTRUCTION
          LM,... LEAKAGE MATRIX FOR RECONSTRUCTION
          PM ... SOURCE VECTOR FOR RECONSTRUCTION
          ] = response_matrix(N, ZON, XDOM, YDOM, ZMAP, QMAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHORTHAND
% XR(m, k, ry, rx)   XL(m, k, ry, rx)   XP(m, ry, rx)
% YR(m, k, ry, rx)   YL(m, k, ry, rx)   YP(m, ry, rx)
% RM(m, k, ry, rx)   LM(m, k, ry, rx)   PM(m, ry, rx)
  
  % INITIALIZE VARIABLES
  NRX = length(XDOM(1, :)); NRY = length(YDOM(1, :));
  [QUAD, chi] = LQN_2D(N);
  [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
  
  % RESPONSE, LEAKAGE MATRICES, AND SOURCE VECTOR INITIALIZATION
  M = N * (N + 2) / 2;
  XR = zeros(M, M, NRY, NRX); YR = zeros(M, M, NRY, NRX);
  XL = zeros(M, M, NRY, NRX); YL = zeros(M, M, NRY, NRX);
  XP = zeros(M, NRY, NRX);    YP = zeros(M, NRY, NRX);
  RM = zeros(M, M, NRY, NRX); LM = zeros(M, M, NRY, NRX);
  PM = zeros(M, NRY, NRX);
  
  % AUXILIARY MATRICES
  XA = zeros(M, M); YA = zeros(M, M); XB = zeros(M, M); YB = zeros(M, M);
  XD = zeros(M, M); YD = zeros(M, M); I = eye(M, M);    H = zeros(M, 1);
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
        
        % XA, YA, XB, YB MATRICES
        if (m <= M / 4) % ROWS
          for k = 1: M  % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            XB(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YB(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > M / 4 && m <= M / 2) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            XB(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            YB(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > M / 2 && m <= 3 * M / 4) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            XB(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YB(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            
          end
        elseif (m > 3 * M / 4 && m <= M) % ROWS
          for k = 1: M % COLUMNS
            XA(m, k) = xvects(m, k, z) * exp(0.5 * st * hx / xvals(k, z));
            XB(m, k) = xvects(m, k, z) * exp(- 0.5 * st * hx / xvals(k, z));
            
            YA(m, k) = yvects(m, k, z) * exp(- 0.5 * st * hy / yvals(k, z));
            YB(m, k) = yvects(m, k, z) * exp(0.5 * st * hy / yvals(k, z));
            
          end
        end
        
        % XD, YD MATRICES
        for k = 1: M % COLUMNS
          miu = QUAD(k, 1); theta = QUAD(k, 2); w = QUAD(k, 3);
          if (k <= M / 4)
            XD(m, k) = 0.25 * c0 * theta * w / (st * leny * (1 - c0));
            if (k == m), XD(m, k) = XD(m,k) + theta / (st * leny); end
         
            YD(m, k) = 0.25 * c0 * miu * w / (st * lenx * (1 - c0));
            if (k == m), YD(m, k) = YD(m,k) + miu / (st * lenx); end
             
          elseif (k > M / 4 && k <= M / 2)
            XD(m, k) = 0.25 * c0 * theta * w / (st * leny * (1 - c0));
            if (k == m), XD(m, k) = XD(m,k) + theta / (st * leny); end
             
            YD(m, k) = - 0.25 * c0 * miu * w / (st * lenx * (1 - c0));
            if (k == m), YD(m, k) = YD(m,k) - miu / (st * lenx); end
              
          elseif (k > M / 2 && k <= 3 * M / 4)
            XD(m, k) = - 0.25 * c0 * theta * w / (st * leny * (1 - c0));
            if (k == m), XD(m, k) = XD(m,k) - theta / (st * leny); end
              
            YD(m, k) = - 0.25 * c0 * miu * w / (st * lenx * (1 - c0));
            if (k == m), YD(m, k) = YD(m,k) - miu / (st * lenx); end
              
          elseif (k > 3 * M / 4 && k <= M)
            XD(m, k) = - 0.25 * c0 * theta * w / (st * leny * (1 - c0));
            if (k == m), XD(m, k) = XD(m,k) - theta / (st * leny); end
              
            YD(m, k) = 0.25 * c0 * miu * w / (st * lenx * (1 - c0));
            if (k == m), YD(m, k) = YD(m,k) + miu / (st * lenx); end
              
          end
          
          XM(m, k) = xvects(m, k, z) * xvals(k, z) * 2 * sinh(0.5 * st * hx / xvals(k, z)) / (st * hx);
          
        end
        
        % H COLUMN VECTOR
        H(m) = Q / (st * (1 - c0));
        
      end
    
      XB_INV = inv(XB); 
      YB_INV = inv(YB);
      
      % CALCULATE THE RESPONSE MATRIX
      XR(:, :, ry, rx) = XA * XB_INV;
      YR(:, :, ry, rx) = YA * YB_INV;
      
      % CALCULATE THE LEAKAGE MATRIX
      XL(:, :, ry, rx) = XD - (XA * (XB_INV * XD));
      YL(:, :, ry, rx) = YD - (YA * (YB_INV * YD));
      
      % CALCULATE THE SOURCE VECTOR
      XP(:, ry, rx) = (I - (XA * XB_INV)) * H;
      YP(:, ry, rx) = (I - (YA * YB_INV)) * H;
      
      
      % CALCULATE THE RESPONSE MATRIX FOR RECONSTRUCTION
      RM(:, :, ry, rx) = XM * XB_INV;
      
      % CALCULATE THE LEAKAGE MATRIX FOR RECONSTRUCTION
      LM(:, :, ry, rx) = XD - XM * (XB_INV * XD);
      
      % CALCULATE THE SOURCE VECTOR FOR RECONSTRUCTION
      PM(:, ry, rx) = (I - XM * XB_INV) * H;
    
    end
  end
  
end