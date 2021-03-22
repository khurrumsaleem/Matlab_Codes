% ANALYTICAL FUNCTION 

function y = myFunc(x, QUAD, c0)
  
  y = 0;
  N = length(QUAD(:,1));
  for m = 1: N
    miu = QUAD(m, 1); w = QUAD(m, 2);
    y = y + 0.5 * c0 * w * x / (miu + x);
  end
  y = y - 1;

end

