% POST - PROCESSING

function POST_2D(N, ZON, XDOM, YDOM, ZMAP, BC, SCALAR_FLUX, X_ANG_FLUX, Y_ANG_FLUX)

  [QUAD, ~] = LQN_2D(N);
  
  % 1. SCALAR FLUX PER REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % CALCULATE THE SCALAR FLUX PER REGIONS
  RX = length(XDOM(1,:)); RY = length(YDOM(1,:));
  mflux_r = zeros(RY, RX);
  for ry = 1: RY
    leny = YDOM(1, ry);
    IB = 0;
    for rx = 1: RX
      lenx = XDOM(1, rx); ncx = XDOM(2, rx); hx = lenx / ncx;
      for i = 1: ncx          
        IB = IB + 1;
        mflux_r(ry, rx) = mflux_r(ry, rx) + leny * hx * SCALAR_FLUX(ry, IB);
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
  
  
  % 2. ABSORPTION RATE PER REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % CALCULATE THE ABSORPTION RATE PER REGION
  r_abs_rate = zeros(RY, RX);
  for ry = 1: RY
    leny = YDOM(1, ry);
    IB = 0;
    for rx = 1: RX
      lenx = XDOM(1, rx); ncx = XDOM(2, rx); z = ZMAP(ry, rx); 
      hx = lenx / ncx; sa = ZON(1, z) - ZON(2, z);
      for i = 1: ncx
        IB = IB + 1;
        r_abs_rate(ry, rx) = r_abs_rate(ry, rx) ...
                             + sa * leny * hx * SCALAR_FLUX(ry, IB);
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
  
  
  % 3. LEAKAGE AT THE BOUNDARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  I = sum(XDOM(2,:)); J = sum(YDOM(2,:));
  M = N * (N + 2) / 2;
  LEAK = zeros(4, 1);
  for ry = 1: RY
    leny = YDOM(1, ry);
    % LEFT
    for m = M / 4 + 1: 3 * M / 4
      miu = QUAD(m, 1); w = QUAD(m, 3);
      LEAK(1) = LEAK(1) + 0.5 * miu * X_ANG_FLUX(m, ry, 1) * leny * w;
    end
    % RIGHT
    for m = 1: M / 4
      miu = QUAD(m, 1); w = QUAD(m, 3);
      LEAK(3) = LEAK(3) + 0.5 * miu * X_ANG_FLUX(m, ry, I + 1) * leny * w;
    end
    for m = 3 * M / 4 + 1: M
      miu = QUAD(m, 1); w = QUAD(m, 3);
      LEAK(3) = LEAK(3) + 0.5 * miu * X_ANG_FLUX(m, ry, I + 1) * leny * w;
    end
  end
  
  for rx = 1: RX
    lenx = XDOM(1, rx);
    % BOTTOM
    for m = M / 2 + 1: M
      theta = QUAD(m, 2); w = QUAD(m, 3);
      LEAK(2) = LEAK(2) + 0.5 * theta * Y_ANG_FLUX(m, 1, rx) * lenx * w;
    end
    % TOP
    for m = 1: M / 2
      theta = QUAD(m, 2); w = QUAD(m, 3);
      LEAK(4) = LEAK(4) + 0.5 * theta * Y_ANG_FLUX(m, J + 1, rx) * lenx * w;
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