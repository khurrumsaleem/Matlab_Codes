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

% HOMOGENEOUS SOLUTION
[xvals, xvects, yvals, yvects] = SPECTRUM_XY(QUAD, chi, ZON);
tol = 10e-10;
for z = 1: NZ
    st = ZON(1,z);
  for k = 1: M
    comp = st .* xvects(:,k,z) ./ xvals(k,z) - A(:,:,z) * xvects(:,k,z);
    for n = 1: M
      assert(comp(n) < tol, 'Homogeneous solution')
    end
  end
end

% PARTICULAR SOLUTION
alfa = zeros(M,M); beta = zeros(M,M); B = zeros(M,1);
st = ZON(1,1); ss = ZON(2,1); c0 = ss/st;
for k = 1: M
  
  % CASE 0
  case0 = 0;
  for m = 1: M
    miu = QUAD(m,1);
    if xvals(k,1) == -miu
        case0 = 1;
    end
  end
  
  if case0 == 0
    % B constant
    num = 0; den = 0;
    for m = 1: M
        miu = QUAD(m,1); w = QUAD(m,3);
        den = den + w*miu*xvects(m,k,1)*xvects(m,k,1);
        num = num + w*xvects(m,k,1)*xvects(m,k,1);
    end
    B(k) = st*num/den;
    
    for m = 1: M
        miu = QUAD(m,1);
        alfa(m,k) = xvects(m,k,1) + 4*xvects(m,k,1)*xvects(m,k,1)/c0 - 4*B(k)*miu*xvects(m,k,1)*xvects(m,k,1)/(c0*st);
        beta(m,k) = B(k) * xvects(m,k,1);
    end
  else
    for m = 1: M
        miu = QUAD(m,1);
        alfa(m,k) = xvects(m,k,1);
        beta(m,k) = st * xvects(m,k,1) / miu;
    end
  end
end

% TEST PARTICULAR SOLUTION
for k = 1: M
  comp0 = st.*alfa(:,k)./xvals(k,1) + beta(:,k) - A(:,:,1)*alfa(:,k) - st.*xvects(:,k,1)./QUAD(:,1);
  comp1 = st.*beta(:,k)./xvals(k,1) - A(:,:,1)*beta(:,k);
  for n = 1: M
    assert(comp0(n) < tol, 'Particular solution');
    assert(comp1(n) < tol, 'Particular solution');
  end
end

