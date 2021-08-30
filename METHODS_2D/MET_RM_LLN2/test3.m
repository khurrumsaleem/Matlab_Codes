N = 2;
ZON = [1; 0.95];
XDOM = [1; 10];
YDOM = [1; 10];
Q = 1;

M = N * (N + 2) / 2;
[QUAD, chi] = LQN_2D(N);
[xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
I = sum(XDOM(2, :)); J = sum(YDOM(2, :));
NZ = length(ZON(1,:));
tol = 10e-12;

st = ZON(1,1); ss = ZON(2,1); c0 = ss/st;
lenx = XDOM(1,1); ncx = XDOM(2,1); hx = lenx/ncx;
leny = YDOM(1,1); ncy = YDOM(2,1); hy = leny/ncy;

% SOLUTION X0
XINV = inv(xvects(:,:,1));
AUX1 = zeros(M,M); AUX2 = zeros(M,M); AUX3 = zeros(M,M); AUX4 = zeros(M,M);
DFLUX0 = rand(M,1); DFLUX1 = rand(M,1); S = zeros(M,1);
for m = 1: M
  for k = 1: M
    k_miu = QUAD(k,1); k_theta = QUAD(k,2);
    AUX1(m,k) = -XINV(m,k)*xvals(m,1)/(st*k_miu);
    AUX2(m,k) = (XINV(m,k)*xvals(m,1)* k_theta / k_miu) / (st * hy);
    AUX3(m,k) = (XINV(m,k)*xvals(m,1)*xvals(m,1)*k_theta/k_miu) * 2 / (st*hx*st*hy);
    AUX4(m,k) = AUX2(m,k);
  end
  S(m) = Q;
end
M1 = xvects(:,:,1)*AUX1*S; M2 = xvects(:,:,1)*AUX2*DFLUX0; 
M3 = xvects(:,:,1)*AUX3*DFLUX1; M4 = xvects(:,:,1)*AUX4*DFLUX1;
PSOL0 = M1 + M2 + M3;
PSOL1 = M4;

% TESTING SOLUTION X0
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2);
  AUX0 = 0; AUX1 = 0;
  for n = 1: M
    nw = QUAD(n,3);
    AUX0 = AUX0 + nw*PSOL0(n);
    AUX1 = AUX1 + nw*PSOL1(n);
  end
  comp0 = 2*miu*PSOL1(m)/(st*hx) + PSOL0(m) - 0.25*c0*AUX0 + theta*DFLUX0(m)/(st*hy) - Q/st;
  comp1 = PSOL1(m) - 0.25*c0*AUX1 + theta*DFLUX1(m)/(st*hy);
  assert(abs(comp0) < tol, 'Error: PSOL0 X0');
  assert(abs(comp1) < tol, 'Error: PSOL1 X0');
end

% SOLUTION X1
AUX1 = zeros(M,M); AUX2 = zeros(M,M); AUX3 = zeros(M,M); AUX4 = zeros(M,M);
AUX5 = zeros(M,M); AUX6 = zeros(M,M); AUX7 = zeros(M,M); AUX8 = zeros(M,M);
DFLUX0 = rand(M,1); DFLUX1 = rand(M,1); S = zeros(M,1);
for m = 1: M
    m_miu = QUAD(m,1); m_theta = QUAD(m,2);
  for k = 1: M
    AUX1(m,k) = 6*m_theta*M1(m,k)/(hy*m_miu); AUX2(m,k) = 6*m_theta*M2(m,k)/(hy*m_miu);
    AUX3(m,k) = 6*m_theta*M3(m,k)/(hy*m_miu); AUX4(m,k) = 6*m_theta*M4(m,k)/(hy*m_miu);
  end
  S(m) = Q;
end
AUX1 = XINV*AUX1; AUX2 = XINV*AUX2; AUX3 = XINV*AUX3; AUX4 = XINV*AUX4;
for m = 1: M
  for k = 1: M
    k_miu = QUAD(k,1); k_theta = QUAD(k,2);
    temp = AUX2(m,k);
    AUX1(m,k) = -AUX1(m,k) * xvals(m,1) / st; AUX2(m,k) = -temp * xvals(m,1) / st;
    AUX3(m,k) = -AUX3(m,k) * xvals(m,1) / st; AUX4(m,k) = -temp * xvals(m,1) / st;
    AUX5(m,k) = -temp * 2 * xvals(m,1) * xvals(m,1) / (st*st*hx);
    AUX6(m,k) = (XINV(m,k)*3*xvals(m,1)*k_theta/k_miu) / (st * hy);
    AUX7(m,k) = (XINV(m,k)*3*xvals(m,1)*k_theta/k_miu) / (st * hy);
    AUX8(m,k) = (XINV(m,k)*3*xvals(m,1)*xvals(m,1)*k_theta/k_miu) * 2 / (st*hx*st*hy);
  end
  S(m) = Q;
end
MM1 = xvects(:,:,1)*AUX1*S; MM2 = xvects(:,:,1)*AUX2*DFLUX0; 
MM3 = xvects(:,:,1)*AUX3*DFLUX1; MM4 = xvects(:,:,1)*AUX4*DFLUX1;
MM5 = xvects(:,:,1)*AUX5*DFLUX1; MM6 = xvects(:,:,1)*AUX6*DFLUX0; 
MM7 = xvects(:,:,1)*AUX7*DFLUX1; MM8 = xvects(:,:,1)*AUX8*DFLUX1;
PSOL0 = MM1 + MM2 + MM3 + MM5 + MM6 + MM8;
PSOL1 = MM4 + MM7;

% TESTING SOLUTION X0
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2);
  AUX0 = 0; AUX1 = 0;
  for n = 1: M
    nw = QUAD(n,3);
    AUX0 = AUX0 + nw*PSOL0(n);
    AUX1 = AUX1 + nw*PSOL1(n);
  end
  comp0 = 2*miu*PSOL1(m)/(st*hx) + PSOL0(m) - 0.25*c0*AUX0 + theta*DFLUX0(m)/(st*hy) - Q/st;
  comp1 = PSOL1(m) - 0.25*c0*AUX1 + theta*DFLUX1(m)/(st*hy);
  assert(abs(comp0) < tol, 'Error: PSOL0 X0');
  assert(abs(comp1) < tol, 'Error: PSOL1 X0');
end


% SOLUTION Y0
YINV = inv(yvects(:,:,1));
AUX1 = zeros(M,M); AUX2 = zeros(M,M); AUX3 = zeros(M,M); AUX4 = zeros(M,M);
DFLUX0 = rand(M,1); DFLUX1 = rand(M,1); S = zeros(M,1);
for m = 1: M
  for k = 1: M
    k_miu = QUAD(k,1); k_theta = QUAD(k,2);
    AUX1(m,k) = -YINV(m,k)*yvals(m,1)/(st*k_theta);
    AUX2(m,k) = (YINV(m,k)*yvals(m,1)*k_miu/k_theta) / (st * hx);
    AUX3(m,k) = (YINV(m,k)*yvals(m,1)*yvals(m,1)*k_miu/k_theta) * 2 / (st*hx*st*hy);
    AUX4(m,k) = AUX2(m,k);
  end
  S(m) = Q;
end
M1 = yvects(:,:,1)*AUX1*S; M2 = yvects(:,:,1)*AUX2*DFLUX0; 
M3 = yvects(:,:,1)*AUX3*DFLUX1; M4 = yvects(:,:,1)*AUX4*DFLUX1;
PSOL0 = M1 + M2 + M3;
PSOL1 = M4;

% TESTING SOLUTION Y0
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2);
  AUX0 = 0; AUX1 = 0;
  for n = 1: M
    nw = QUAD(n,3);
    AUX0 = AUX0 + nw*PSOL0(n);
    AUX1 = AUX1 + nw*PSOL1(n);
  end
  comp0 = 2*theta*PSOL1(m)/(st*hy) + PSOL0(m) - 0.25*c0*AUX0 + miu*DFLUX0(m)/(st*hx) - Q/st;
  comp1 = PSOL1(m) - 0.25*c0*AUX1 + miu*DFLUX1(m)/(st*hx);
  assert(abs(comp0) < tol, 'Error: PSOL0 X0');
  assert(abs(comp1) < tol, 'Error: PSOL1 X0');
end
