% CALCULATE THE ROOT BETWEEN a AND b USING BISECT METHOD

function root = myRootFunc(a, b, QUAD, c0, tol)

  yb = myFunc(b, QUAD, c0);
  while (abs(b - a) > tol)
      c = (a + b) / 2; yc = myFunc(c, QUAD, c0);
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