% TESTING QUADRATURE PROPERTIES
N = 8;
M = N * (N + 2) / 2;
[QUAD, chi] = LQN_2D(N);

% SUM(w * miu)
sum = 0;
for m = 1: M
  miu = QUAD(m, 1); w = QUAD(m, 3);
  sum = sum + w * miu;
end
fprintf("SUM(w * miu) = %.5e for N = %d\n", sum, N);

% SUM(w * theta)
sum = 0;
for m = 1: M
  theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * theta;
end
fprintf("SUM(w * theta) = %.5e for N = %d\n", sum, N);

% SUM(w * theta^2)
sum = 0;
for m = 1: M
  theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * theta * theta;
end
fprintf("SUM(w * theta^2) = %.5f for N = %d\n", sum, N);

% SUM(w * miu * theta)
sum = 0;
for m = 1: M
  miu = QUAD(m, 1); theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * miu * theta;
end
fprintf("SUM(w * miu * theta) = %.5e for N = %d\n", sum, N);

% SUM(w * miu * theta^2)
sum = 0;
for m = 1: M
  miu = QUAD(m, 1); theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * miu * theta * theta;
end
fprintf("SUM(w * miu * theta^2) = %.5f for N = %d\n", sum, N);

% SUM(w * miu * xvects * xvects)
[xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);

for k = 1: M
  sum = 0;
  for m = 1: M
    miu = QUAD(m, 1); theta = QUAD(m, 2); w = QUAD(m, 3);
    sum = sum + w * miu * xvects(m,k,1) * xvects(m,k,1);    
  end
  fprintf("SUM(w * miu * xvect^2) = %.5f for N = %d\n", sum, N);
end

for k = 1: M
  sum = 0;
  for m = 1: M
    miu = QUAD(m, 1); theta = QUAD(m, 2); w = QUAD(m, 3);
    sum = sum + w * theta * xvects(m,k,1) / miu;    
  end
  fprintf("miu = %.5f, val = %.5f\n", miu, xvals(k,1));
  fprintf("SUM(w * theta * xvect / miu) = %.5f for N = %d\n", sum, N);
end
