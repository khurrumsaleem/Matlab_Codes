% CONSTRUCTION OF THE RESPONSE MATRIX R AND SOURCE MATRIX P

function [R, P] = response_matrix(QUAD, ZON, DOM)
% R(m, k, nr)   % P(m, nr)
  
  % INITIALIZE VARIABLES
  NR = length(DOM(1, :)); N = length(QUAD(:,1));
  [vals, vects] = spectrum(QUAD, ZON);
  R = zeros(N, N, NR); P = zeros(N, NR);
  A = zeros(N, N); B = zeros(N, N); I = eye(N, N); S = zeros(N, 1);
  
  % DEFINE EACH MATRIX PER REGION
  for nr = 1: NR
    
    % CONSTRUCT AUXILIARY MATRICES A AND B
    len = DOM(1, nr); ntc = DOM(2, nr); h = len / ntc;
    z = DOM(3, nr); st = ZON(1, z); ss = ZON(2, z); c0 = ss/st;
    Q = DOM(4, nr);
    for m = 1: N
      if m <= N/2
        for k = 1: N
          A(m, k) = vects(m, k) * exp(st * h / vals(k));
          B(m, k) = vects(m, k);
        end
      else
        for k = 1: N
          A(m, k) = vects(m, k);
          B(m, k) = vects(m, k) * exp(st * h / vals(k));
        end
      end
      S(m) = Q / (st * (1 - c0));
    end
    
    % OBTAIN R AND P MATRICES
    AUX = A * inv(B);
    R(:, :, nr) = AUX;
    P(:, nr) = (I - AUX) * S;
    
  end
  
end