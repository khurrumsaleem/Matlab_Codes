% TESTING QUADRATURE PROPERTIES
N = 4;
M = N * (N + 2) / 2;
[QUAD, chi] = LQN_2D(N);

% SUM(w * miu)
sum = 0;
for m = 1: M
  miu = QUAD(m, 1); w = QUAD(m, 3);
  sum = sum + w * miu;
end
fprintf("SUM(w * miu) = %.5f for N = %d\n", sum, N);

% SUM(w * theta)
sum = 0;
for m = 1: M
  theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * theta;
end
fprintf("SUM(w * theta) = %.5f for N = %d\n", sum, N);

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
fprintf("SUM(w * miu * theta) = %.5f for N = %d\n", sum, N);

% SUM(w * miu * theta^2)
sum = 0;
for m = 1: M
  miu = QUAD(m, 1); theta = QUAD(m, 2); w = QUAD(m, 3);
  sum = sum + w * miu * theta * theta;
end
fprintf("SUM(w * miu * theta^2) = %.5f for N = %d\n", sum, N);