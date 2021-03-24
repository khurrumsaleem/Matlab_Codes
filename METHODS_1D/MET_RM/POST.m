% POST-PROCESSING

% OUTPUT: AV_SCALAR_FLUX = AVERAGE SCALAR FLUX
%         AV_ANGULAR_FLUX = AVERAGE ANGULAR FLUX
%         ABS_R = ABSORPTION RATE PER REGION
%         LEAK = LEAKAGES AT THE BOUNDARIES
%         X = SPATIAL DISCRETIZATION

function [AV_SCALAR_FLUX, AV_ANGULAR_FLUX, ABS_R, LEAK, X] = ...
         POST(N, ZON, DOM, BC, SCALAR_FLUX, ANGULAR_FLUX)

  % VARIABLES
  ntc = sum(DOM(2, :));   QUAD = GAUSS_QUAD(N);
  NR = length(DOM(1, :)); X = zeros(ntc + 1, 1);
  [vals, vects] = spectrum(QUAD, ZON);
  R = zeros(N, N, NR); P = zeros(N, NR);
  C = zeros(N, N); B = zeros(N, N); I = eye(N, N); S = zeros(N, 1);
  
  % PLOT SCALAR FLUX
  IB = 0;
  for r = 1: NR
    nc = DOM(2,r);
    len = DOM(1, r); h = len / nc;
    for i = 1: nc
      IB = IB + 1; IF = IB + 1; 
      X(IF) = X(IB) + h;
    end
  end
  plot(X, SCALAR_FLUX);
  
  for nr = 1: NR
    % CONSTRUCT AUXILIARY MATRICES C AND B
    len = DOM(1, nr); ntc = DOM(2, nr); h = len / ntc;
    z = DOM(3, nr); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
    Q = DOM(4, nr);
    for m = 1: N
      if m <= N/2
        for k = 1: N
          B(m, k) = vects(m, k);
          C(m, k) = vects(m, k) * vals(k) * (exp(st * h / vals(k)) - 1) / (st * h);
        end
      else
        for k = 1: N
          B(m, k) = vects(m, k) * exp(st * h / vals(k));
          C(m, k) = vects(m, k) * vals(k) * (exp(st * h / vals(k)) - 1) / (st * h);
        end
      end
      S(m) = Q / (st * (1 - c0));
    end
    
    % OBTAIN R AND P MATRICES
    AUX = C * inv(B);
    R(:, :, nr) = AUX;
    P(:, nr) = (I - AUX) * S;
  end
  
  % AVERAGE ANGULAR FLUX CALCULATION
  AV_ANGULAR_FLUX = zeros(ntc, N);
  AV_SCALAR_FLUX = zeros(ntc, 1);
  IB = 0;
  for r = 1: NR
    nc = DOM(2, r);
    RM = R(:, :, r); PM = P(:, r);
    for i = 1: nc
      IB = IB + 1; IF = IB + 1;
      IN = [ANGULAR_FLUX(IB, 1: N / 2) ANGULAR_FLUX(IF, N / 2 + 1: N)]';
      AV_ANGULAR_FLUX(IB, :) = RM * IN + PM;
      aux = 0;
      for m = 1: N
        w = QUAD(m, 2);
        aux = aux + AV_ANGULAR_FLUX(IB, m) * w;
      end
      AV_SCALAR_FLUX(IB) = 0.5 * aux;
    end
  end
  
  % PRINT ABSORPTION RATE PER REGION
  fprintf("\nR\t\tABSORPTION_RATE\n");
  ABS_R = zeros(NR, 1);
  IB = 0;
  for r = 1: NR
    nc = DOM(2, r); z = DOM(3, r); st = ZON(1, z); ss = ZON(2, z);
    sa = st - ss; len = DOM(1, r); h = len / nc;
    for i = 1: nc
      IB = IB + 1; 
      for m = 1: N
        ABS_R(r) = ABS_R(r) + sa * h * AV_SCALAR_FLUX(IB);
      end
    end
    fprintf("%d\t\t%.4e\n", r, ABS_R(r));
  end
  
  % PRINT BONDARY LEAKAGES
  fprintf("\n\t\tLEAKAGES\n");
  LEFT = 0.0; RIGHT = 0.0;
  for m = 1: N / 2
    w = QUAD(m, 2); miu = QUAD(m, 1);
    RIGHT = RIGHT + 0.5 * h * w * miu * ANGULAR_FLUX(ntc + 1, m);
    w = QUAD(N / 2 + m, 2); miu = QUAD(N / 2 + m, 1);
    LEFT = LEFT + 0.5 * h * w * miu * ANGULAR_FLUX(1, N / 2 + m);
  end
  if (BC(1) == -1)
    LEFT = 0;
    fprintf("LEFT\t-\n");
  else
    fprintf("LEFT\t%.4e\n", LEFT);
  end
  if (BC(2) == -1)
    RIGHT = 0;
    fprintf("RIGHT\t-\n");
  else
    fprintf("RIGHT\t%.4e\n", RIGHT);
  end
  LEAK = [LEFT, RIGHT];
  
end