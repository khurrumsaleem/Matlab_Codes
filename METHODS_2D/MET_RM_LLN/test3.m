[QUAD, chi] = LQN_2D(N);

M = N * (N + 2) / 2; I = sum(XDOM(2, :)); J = sum(YDOM(2, :));
NZ = length(ZON(1,:));
tol = 10e-10;

A = zeros(M, M, NZ);
for z = 1: NZ
  st = ZON(1,z); ss = ZON (2,z); c0 = ss/st;
  for m = 1: M
    miu = QUAD(m,1);
    for n = 1: M
      w = QUAD(n,3);
      A(m,n) =0.25*st*c0*w/miu;
      if (m == n)
        A(m,n) = 0.25*st*c0*w/miu - st/miu ;
      end
    end
  end
end

% PARTICULAR SOLUTION X0
DFLUX0 = rand(M,1); DFLUX1 = rand(M,1);
PSOL0 = zeros(M,1); PSOL1 = zeros(M,1);
st = ZON(1,1); ss = ZON(2,1); c0 = ss/st; hi = rand(1); hj = rand(1);
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2); w = QUAD(m,3);
  AUX0 = 0; AUX1 = 0; AUX2 = 0; AUX4 = 0;
  for k = 1: M
    kmiu = QUAD(k,1); ktheta = QUAD(k,2); kw = QUAD(k,3);
    AUX0 = AUX0 + 0.25*c0*kw*ktheta*DFLUX0(k)/(st*hj*(1-c0));
    AUX1 = AUX1 + 0.5*c0*miu*kw*ktheta*DFLUX1(k)/(st*hi*st*hj*(1-c0));
    AUX2 = AUX2 + 0.5*c0*kw*kmiu*ktheta*DFLUX1(k)/(st*hi*st*hj*(1-c0));
    AUX4 = AUX4 + 0.25*c0*kw*ktheta*DFLUX1(k)/(st*hj*(1-c0));
  end
  PSOL0(m) = -theta*DFLUX0(m)/(st*hj) - AUX0 + 2*miu*theta*DFLUX1(m)/(st*hi*st*hj) ...
             + AUX1 + AUX2;
  PSOL1(m) = -theta*DFLUX1(m)/(st*hj) - AUX4;
end
% TESTING PARTICULAR SOLUTION Y0
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2);
  AUX0 = 0; AUX1 = 0;
  for n = 1: M
    nw = QUAD(n,3);
    AUX0 = AUX0 + nw*PSOL0(n);
    AUX1 = AUX1 + nw*PSOL1(n);
  end
  comp0 = 2*miu*PSOL1(m)/(st*hi) + PSOL0(m) - 0.25*c0*AUX0 + theta*DFLUX0(m)/(st*hj);
  comp1 = PSOL1(m) - 0.25*c0*AUX1 + theta*DFLUX1(m)/(st*hj);
  assert(comp0 < tol, 'Error: PSOL0 X0');
  assert(comp1 < tol, 'Error: PSOL1 X0');
end


% PARTICULAR SOLUTION Y0
DFLUX0 = rand(M,1); DFLUX1 = rand(M,1);
PSOL0 = zeros(M,1); PSOL1 = zeros(M,1);
st = ZON(1,1); ss = ZON(2,1); c0 = ss/st; hi = rand(1); hj = rand(1);
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2); w = QUAD(m,3);
  AUX0 = 0; AUX1 = 0; AUX2 = 0; AUX4 = 0;
  for k = 1: M
    kmiu = QUAD(k,1); ktheta = QUAD(k,2); kw = QUAD(k,3);
    AUX0 = AUX0 + 0.25*c0*kw*kmiu*DFLUX0(k)/(st*hi*(1-c0));
    AUX1 = AUX1 + 0.5*c0*theta*kw*kmiu*DFLUX1(k)/(st*hi*st*hj*(1-c0));
    AUX2 = AUX2 + 0.5*c0*kw*kmiu*ktheta*DFLUX1(k)/(st*hi*st*hj*(1-c0));
    AUX4 = AUX4 + 0.25*c0*kw*kmiu*DFLUX1(k)/(st*hi*(1-c0));
  end
  PSOL0(m) = -miu*DFLUX0(m)/(st*hi) - AUX0 + 2*miu*theta*DFLUX1(m)/(st*hi*st*hj) ...
             + AUX1 + AUX2;
  PSOL1(m) = -miu*DFLUX1(m)/(st*hi) - AUX4;
end
% TESTING PARTICULAR SOLUTION Y0
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(m,2);
  AUX0 = 0; AUX1 = 0;
  for n = 1: M
    nw = QUAD(n,3);
    AUX0 = AUX0 + nw*PSOL0(n);
    AUX1 = AUX1 + nw*PSOL1(n);
  end
  comp0 = 2*theta*PSOL1(m)/(st*hj) + PSOL0(m) - 0.25*c0*AUX0 + miu*DFLUX0(m)/(st*hi);
  comp1 = PSOL1(m) - 0.25*c0*AUX1 + miu*DFLUX1(m)/(st*hi);
  assert(comp0 < tol, 'Error: PSOL0 Y0');
  assert(comp1 < tol, 'Error: PSOL1 Y0');
end


% PARTICULAR SOLUTION
[xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
INV = inv(xvects);

st = 1; hj = 1;
NABLA = zeros(M,M); F0 = zeros(M,M); F1 = zeros(M,M);
for m = 1: M
  miu = QUAD(m,1); theta = QUAD(M,2);
  for k = 1: M
    NABLA(m,k) = 6*theta*xvects(m,k)/(hj*miu);
    if (xvals(m,1) == xvals(k,1))
      F0(m,k) = 0;
      F1(m,k) = 1;
    else
      F0(m,k) = xvals(m,1)*xvals(k,1)/(st*(xvals(m,1) - xvals(k,1)));
      F1(m,k) = 0;
    end
  end
end
B = INV*NABLA;
landa = xvects*(B.*F0);
gama = xvects*(B.*F1);

% TEST PARTICULAR SOLUTION
for k = 1: M
  comp0 = st.*landa(:,k)./xvals(k,1) + gama(:,k) - A(:,:,1)*landa(:,k) - NABLA(:,k);
  comp1 = st.*gama(:,k)./(xvals(k,1)) - A(:,:,1)*gama(:,k);
  for n = 1: M
    assert(comp0(n) < tol, 'Particular solution');
    assert(comp1(n) < tol, 'Particular solution');
  end
end