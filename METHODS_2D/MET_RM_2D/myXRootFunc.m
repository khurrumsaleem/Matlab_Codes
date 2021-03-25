% CALCULATE THE ROOT BETWEEN a AND b FOR THE DISPERSION LAW IN X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function root = myXRootFunc(a, b, QUAD, c0)
 
  yb = myXFunc(b, QUAD, c0);
  while (abs(b - a) > 1e-12)
      c = (a + b) / 2; yc = myXFunc(c, QUAD, c0);
      if yc == 0
          a = c; b = c;
      elseif yb * yc > 0
          b = c; yb = yc;
      else
          a = c;
      end
  end
  root = (a + b) / 2;

end