% DISPERSION LAW IN Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = myYFunc(x, QUAD, c0)

  y = 0;
  M = length(QUAD(:, 1));
  for m = 1: M
    theta = QUAD(m, 2); w = QUAD(m, 3);
    y = y + w / (theta + x);
  end
  y = 0.25 * c0 * x * y - 1;

end

