% DISPERSION LAW IN X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = myXFunc(x, QUAD, c0)

  y = 0; M = length(QUAD(:, 1));
  for m = 1: M
    miu = QUAD(m, 1); w = QUAD(m, 3);
    y = y + w / (miu + x);
  end
  y = 0.25 * c0 * x * y - 1;

end

