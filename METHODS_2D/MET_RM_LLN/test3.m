[QUAD, chi] = LQN_2D(N);

M = N * (N + 2) / 2; I = sum(XDOM(2, :)); J = sum(YDOM(2, :));
NZ = length(ZON(1,:));

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

% PARTICULAR SOLUTION
[xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
INV = inv(xvects);
tol = 10e-10;
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