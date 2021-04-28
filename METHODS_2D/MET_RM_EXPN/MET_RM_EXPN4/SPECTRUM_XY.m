% EIGENVALUES AND EIGENVECTORS CALCULATION FOR THE TRANSVERSAL INTEGRATED
% SN EQUATIONS FOR PROBLEMS IN 2D

% OUTPUT: xvals = EIGENVALUES FOR OPERATOR IN X
%         xvects = EIGENVECTORS FOR OPERATOR IN X
%         yvals = EIGENVALUES FOR OPERATOR IN Y
%         yvects = EIGENVECTORS FOR OPERATOR IN Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON)

% SHORTHANDS
% xvals(k, z)     yvals(k, z)
% xvects(m, k, z) yvects(m, k, z)
  
  % VARIABLES
  M = length(QUAD(:, 1)); N = sqrt(2 * M + 1) - 1; NZ = length(ZON(1, :));
  xvals = zeros(M, NZ); xvects = zeros(M, M, NZ);
  yvals = zeros(M, NZ); yvects = zeros(M, M, NZ);
  
  for z = 1: NZ
    st = ZON(1, z); ss = ZON(2, z); c0 = ss / st;
    if (c0 ~= 0)
        
      % EIGENVALUES CALCULATION FOR X
      xvals1 = zeros(N, 1); xvals2 = zeros(M - N, 1);
      
      % BY DISPERSION LAW
      for i = 1: N / 2
        chi_i = chi(i);
        if i == N / 2, chi_f = 10; else, chi_f = chi(i + 1); end
        h = (chi_f - chi_i) / (10 ^ N); a = chi_i + h; b = chi_f - h;
        r = myXRootFunc(a, b, QUAD, c0); xvals1(i) = r;
      end
      % BY ZERO NORMALIZATION
      k = 0;
      for i = 1: N / 2
        val = chi(i); xmult = 0;
        for m = 1: M
          miu = QUAD(m, 1);
          if ( val == miu)
            xmult = xmult + 1;
          end
        end
        xmult = xmult - 1;
        while (xmult > 0)
          k = k + 1; xvals2(k) = val; xmult = xmult - 1;
        end
      end
      
      % EIGENVALUES CALCULATION FOR Y
      yvals1 = zeros(N,1); yvals2 = zeros(M-N,1);
      
      % BY DISPERSION LAW
      for i = 1: N / 2
        chi_i = chi(i);
        if i == N / 2, chi_f = 10; else, chi_f = chi(i + 1); end
        h = (chi_f - chi_i) / (10 ^ N); a = chi_i + h; b = chi_f - h;
        r = myYRootFunc(a, b, QUAD, c0); yvals1(i) = r;
      end
      % BY ZERO NORMALIZATION
      k = 0;
      for i = 1: N / 2
        val = chi(i); ymult = 0;
        for m = 1: M
          theta = QUAD(m, 1);
          if ( val == theta)
            ymult = ymult + 1;
          end
        end
        ymult = ymult - 1;
        while (ymult > 0)
          k = k + 1; yvals2(k) = val; ymult = ymult - 1;
        end
      end
      
      % ORDERING THE EIGENVALUES OF X
      if (c0 ~= 0), xvals1(:) = sort(xvals1(:), 'descend'); end
      for i = 1: (M - N) / 2
        if (i <= N / 2), xvals1(N / 2 + i) = - xvals1(i); end
        xvals2((M - N) / 2 + i) = - xvals2(i);
      end
      xvals(:, z) = [xvals2(:); xvals1(:)];
      
      % % ORDERING THE EIGENVALUES OF Y
      if (c0 ~= 0), yvals1(:) = sort(yvals1(:), 'descend'); end
      for i = 1: (M - N) / 2
        if (i <= N / 2), yvals1(N / 2 + i) = - yvals1(i); end
        yvals2((M - N) / 2 + i) = - yvals2(i);
      end
      yvals(:, z) = [yvals2(:); yvals1(:)];
      
    else % FOR C0 = 0
        
      % EIGENVALUES FOR C0 = 0
      for i = 1: M
        miu = QUAD(i, 1); theta = QUAD(i, 2);
        xvals(i, z) = - miu;
        yvals(i, z) = - theta;
      end
      
    end
    
    if (c0 ~= 0)
        
      % EIGENVECTOR CALCULATIONS BY DISPERSION LAW IN X
      xvects1 = zeros(M, N);
      for i = 1: N
        for m = 1: M
          miu = QUAD(m, 1);
          xvects1(m, i) = 0.25 * c0 * xvals1(i) / (miu + xvals1(i));
        end
      end
      % BY ZERO NORMALIZATION
      xvects2 = zeros(M, M - N); aux = zeros(M, 1);
      for i = 1: M - N
        val = xvals2(i);
        for m = 1: M
          mw = QUAD(m, 3);
          mmiu = QUAD(m, 1);
          if (val == - mmiu && aux(m) == 0)
            aux(m) = 1;
            for n = 1: M
              nw = QUAD(n, 3);
              nmiu = QUAD(n, 1);
              if(val == - nmiu && aux(n) == 0)
                aux(n) = 1;
                xvects2(n, i) = - mw / nw;
                xvects2(m, i) = 1;
                aux(m) = 0;
                break;
              end
            end
            break;
          end
        end
      end
      xvects(:, :, z) = [xvects2(:, :), xvects1(:, :)];
      
      % EIGENVECTOR CALCULATIONS BY DISPERSION LAW IN Y
      yvects1 = zeros(M, N);
      for i = 1: N
        for m = 1: M
          theta = QUAD(m, 2);
          yvects1(m, i) = 0.25 * c0 * yvals1(i) / (theta + yvals1(i));
        end
      end
      % BY ZERO NORMALIZATION
      yvects2 = zeros(M, M - N); aux = zeros(M, 1);
      for i = 1: M - N
        val = yvals2(i);
        for m = 1: M
          mw = QUAD(m, 3);
          mtheta = QUAD(m, 2);
          if (val == - mtheta && aux(m) == 0)
            aux(m) = 1;
            for n = 1: M
              nw = QUAD(n, 3);
              ntheta = QUAD(n, 2);
              if(val == - ntheta && aux(n) == 0)
                aux(n) = 1;
                yvects2(n, i) = - mw / nw;
                yvects2(m, i) = 1;
                aux(m) = 0;
                break;
              end
            end
            break;
          end
        end
      end
      yvects(:, :, z) = [yvects2(:, :), yvects1(:, :)];
      
    else % FOR C0 = 0
        
      % EIGENVECTORS FOR C0 = 0
      for i = 1: M
        for m = 1: M
          if (m == i)
            xvects(m, i, z) = 1.0; yvects(m, i, z) = 1.0;
          else
            xvects(m, i, z) = 0.0; yvects(m, i, z) = 0.0;
          end
        end
      end
      
    end
  end

end