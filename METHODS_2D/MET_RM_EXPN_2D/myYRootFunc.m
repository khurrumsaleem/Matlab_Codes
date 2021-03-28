% CALCULATE THE ROOT BETWEEN a AND b FOR THE DISPERSION LAW IN X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = myYRootFunc(a, b, QUAD, c0)
 
  yb = myYFunc(b,QUAD,c0);
  while (abs(b - a) > 1e-12)
      c = (a + b) / 2; yc = myYFunc(c, QUAD, c0);
      if yc == 0
          a = c; b = c;
      elseif yb * yc > 0
          b = c; yb = yc;
      else
          a = c;
      end
  end
  r = (a + b) / 2;

end