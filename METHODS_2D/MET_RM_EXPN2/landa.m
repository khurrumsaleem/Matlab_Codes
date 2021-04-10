% CALCULATE LANDA FOR EACH NODE

function [LANDA0, ... % LANDA IN THE CURRENT NODE
          LANDAL, ... % LANDA IN THE LEFT NODE
          LANDAR, ... % LANDA IN THE RIGHT NODE
          LANDAB, ... % LANDA IN THE BOTTOM NODE
          LANDAT  ... % LANDA IN THE TOP NODE
          ] = landa(XDOM, YDOM, ZON, ZMAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHORTHAND
% LANDA0(JB, IB)
% LANDAL(JB, IB)     LANDAR(JB, IB)
% LANDAB(JB, IB)     LANDAT(JB, IB)

  RX = length(XDOM(1, :)); RY = length(YDOM(1,:));
  I = sum(XDOM(2, :));     J = sum(YDOM(2, :));
  
  % LANDA IN THE CURRENT NODE
  LANDA0 = zeros(J, I);
  JB = 0;
  for ry = 1: RY
    ncy = YDOM(2, ry);
    for j = 1: ncy
      JB = JB + 1;
      IB = 0;
      for rx = 1: RX
        z = ZMAP(ry, rx); ss = ZON(2,z);
        ncx = XDOM(2, rx);
        for i = 1: ncx
          IB = IB + 1;
          LANDA0(JB, IB) = ss;
        end
      end
    end
  end
  
  LANDAL = zeros(J, I); LANDAR = zeros(J, I);
  LANDAB = zeros(J, I); LANDAT = zeros(J, I);
  JB = 0;
  for ry = 1: RY
    ncy = YDOM(2, ry);
    for j = 1: ncy
      JB = JB + 1;
      IB = 0;
      for rx = 1: RX
        ncx = XDOM(2, rx);
        for i = 1: ncx
          IB = IB + 1;
          
          % LANDA IN THE LEFT NODE
          if (rx == 1 && i == 1)
            LANDAL(JB, IB) = LANDA0(JB, IB);
          else
            LANDAL(JB, IB) = LANDA0(JB, IB - 1);
          end
          
          % LANDA IN THE RIGHT NODE
          if (rx == RX && i == ncx)
            LANDAR(JB, IB) = LANDA0(JB, IB);
          else
            LANDAR(JB, IB) = LANDA0(JB, IB + 1);
          end
          
          % LANDA IN THE BOTTOM NODE
          if (ry == 1 && j == 1)
            LANDAB(JB, IB) = LANDA0(JB, IB);
          else
            LANDAB(JB, IB) = LANDA0(JB - 1, IB);
          end
          
          % LANDA IN THE TOP NODE
          if (ry == RY && j == ncy)
            LANDAT(JB, IB) = LANDA0(JB, IB);
          else
            LANDAT(JB, IB) = LANDA0(JB + 1, IB);
          end
          
        end
      end
    end
  end

end