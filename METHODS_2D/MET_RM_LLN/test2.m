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
          
        
      