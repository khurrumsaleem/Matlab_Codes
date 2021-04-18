% POST - PROCESSING. VERSION 4

function POST_2D(N, ZON, XDOM, YDOM, ZMAP, BC, SCALAR_FLUX, X_ANG_FLUX, Y_ANG_FLUX)

  [QUAD, ~] = LQN_2D(N);
  
  % 1. PLOT SCALAR FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    I = sum(XDOM(2, :)); J = sum(YDOM(2, :));
    RX = length(XDOM(1, :)); RY = length(YDOM(1, :));
    if (J > 1 && I > 1)
      % CONSTRUCTION OF X AXIS
      X = zeros(I + 1, 1); IB = 0;
      for rx = 1: RX
        lenx = XDOM(1, rx); ncx = XDOM(2, rx); xh = lenx / ncx;
        for i = 1: ncx
          IB = IB + 1; IF = IB + 1;
          X(IF, 1) = X(IB, 1) + xh;    
        end
      end
  
      % CONSTRUCTION OF Y AXIS
      Y = zeros(J + 1, 1); JB = 0;
      for ry = 1: RY
        leny = YDOM(1, ry); ncy = YDOM(2, ry); yh = leny / ncy;
        for i = 1: ncy
          JB = JB + 1; JF = JB + 1;
          Y(JF, 1) = Y(JB, 1) + yh;    
        end
      end
    
      % PLOT SCALAR FLUX
      xm = zeros(length(X) - 1, 1); ym = zeros(length(Y) - 1, 1);
      for i = 1: length(X) - 1
        xm(i) = 0.5 * (X(i + 1) + X(i));
      end
      for j = 1: length(Y) - 1
        ym(j) = 0.5 * (Y(j + 1) + Y(j));
      end
      [x, y] = meshgrid(xm, ym); 
      surf(x, y, SCALAR_FLUX); 
      hold on; 
      contourf(x, y, SCALAR_FLUX);
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % 2. SCALAR FLUX PER REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % CALCULATE THE SCALAR FLUX PER REGIONS
  mflux_r = zeros(RY, RX);
  JB = 0;
  for ry = 1: RY
    leny = YDOM(1, ry); ncy = YDOM(2, ry); hy = leny / ncy;
    for j = 1: ncy
      JB = JB + 1;
      IB = 0;
      for rx = 1: RX
        lenx = XDOM(1, rx); ncx = XDOM(2, rx); hx = lenx / ncx;
        for i = 1: ncx
          IB = IB + 1;
          mflux_r(ry, rx) = mflux_r(ry, rx) + hy * hx * SCALAR_FLUX(JB, IB);
        end
      end
    end
  end
  
  % PRINT SCALAR FLUX PER REGIONS
  fprintf("\nSCALAR FLUX PER REGION:\n");
  for rx = 1: RX
     if(rx == 1), fprintf("RY RX"); end 
     fprintf("\t\t%d\t\t", rx);
  end
  fprintf("\n");
  for ry = 1: RY
    fprintf("%d", ry);
    leny = YDOM(1, ry);
    for rx = 1: RX
      lenx = XDOM(1,rx); area = leny * lenx;
      mflux_r(ry, rx) = mflux_r(ry, rx) / area;
      fprintf("\t\t%.5E", mflux_r(ry, rx));
    end
    fprintf("\n");
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % 3. ABSORPTION RATE PER REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % CALCULATE THE ABSORPTION RATE PER REGION
  r_abs_rate = zeros(RY, RX); JB = 0;
  for ry = 1: RY
    leny = YDOM(1, ry); ncy = YDOM(2, ry); hy = leny / ncy;
    for j = 1: ncy
      JB = JB + 1;
      IB = 0;
      for rx = 1: RX
        lenx = XDOM(1, rx); ncx = XDOM(2, rx); z = ZMAP(ry, rx); 
        hx = lenx / ncx; sa = ZON(1, z) - ZON(2, z);
        for i = 1: ncx
          IB = IB + 1;
          r_abs_rate(ry, rx) = r_abs_rate(ry, rx) ...
                               + sa * hy * hx * SCALAR_FLUX(JB, IB);
        end
      end
    end
  end
  
  % PRINT ABSORPTION RATE
  fprintf("\nABSORPTION RATE PER REGION:\n");
  for rx = 1: RX
     if(rx == 1), fprintf("RY RX"); end 
     fprintf("\t\t%d\t\t", rx);
  end
  fprintf("\n");
  for ry = 1: RY
    fprintf("%d", ry);
    for rx = 1: RX
      fprintf("\t\t%.5E", r_abs_rate(ry, rx));
    end
    fprintf("\n");
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % 4. LEAKAGE AT THE BOUNDARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  M = N * (N + 2) / 2;
  LEAK = zeros(4, 1); JB = 0;
  for ry = 1: RY
    leny = YDOM(1, ry); ncy = YDOM(2, ry); hy = leny / ncy;
    for j = 1: ncy
      JB = JB + 1;
      % LEFT
      for m = M / 4 + 1: 3 * M / 4
        miu = QUAD(m, 1); w = QUAD(m, 3);
        LEAK(1) = LEAK(1) + 0.5 * miu * X_ANG_FLUX(m, JB, 1) * hy * w;
      end
      % RIGHT
      for m = 1: M / 4
        miu = QUAD(m, 1); w = QUAD(m, 3);
        LEAK(3) = LEAK(3) + 0.5 * miu * X_ANG_FLUX(m, JB, I + 1) * hy * w;
      end
      for m = 3 * M / 4 + 1: M
        miu = QUAD(m, 1); w = QUAD(m, 3);
        LEAK(3) = LEAK(3) + 0.5 * miu * X_ANG_FLUX(m, JB, I + 1) * hy * w;
      end
    end
  end
  IB = 0;
  for rx = 1: RX
    lenx = XDOM(1, rx); ncx = XDOM(2, rx); hx = lenx / ncx;
    for i = 1: ncx
      IB = IB + 1;
      % BOTTOM
      for m = M / 2 + 1: M
        theta = QUAD(m, 2); w = QUAD(m, 3);
        LEAK(2) = LEAK(2) + 0.5 * theta * Y_ANG_FLUX(m, 1, IB) * hx * w;
      end
      % TOP
      for m = 1: M / 2
        theta = QUAD(m, 2); w = QUAD(m, 3);
        LEAK(4) = LEAK(4) + 0.5 * theta * Y_ANG_FLUX(m, J + 1, IB) * hx * w;
      end
    end
  end

  fprintf("\nLEAKAGE:\n");
  fprintf("LEFT\t\t\tBOTTOM\t\t\tRIGHT\t\t\tTOP\n");
  for i = 1: 4
     if (BC(i) == -1), fprintf("-\t\t\t\t");
     else, fprintf("% .4E\t\t", LEAK(i)); end
  end
  fprintf("\n");
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end