% CALCULATE LAMDA FOR EACH REGION

function [LAMDA0, ... % LAMDA IN THE CURRENT REGION
          LAMDAL, ... % LAMDA IN THE LEFT REGION
          LAMDAR, ... % LAMDA IN THE RIGHT REGION
          LAMDAB, ... % LAMDA IN THE BOTTOM REGION
          LAMDAT  ... % LAMDA IN THE TOP REGION
          ] = lamda(XDOM, YDOM, ZON, ZMAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHORTHAND
% LAMDA0(ry,rx)
% LAMDAL(ry,rx)     LAMDAR(ry,rx)
% LAMDAB(ry,rx)     LAMDAT(ry,rx)

  RX = length(XDOM(1, :)); RY = length(YDOM(1,:));
  
  % LAMDA IN THE CURRENT NODE
  LAMDA0 = zeros(RY, RX);
  for ry = 1: RY
    for rx = 1: RX
      z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z);
      LAMDA0(ry, rx) = st - ss;
    end
  end
  
  % LAMDA IN THE LEFT REGION
  LAMDAL = zeros(RY, RX);
  for ry = 1: RY
    for rx = 1: RX
      if (rx == 1)
        LAMDAL(ry, rx) = LAMDA0(ry, rx);
      else
        LAMDAL(ry, rx) = LAMDA0(ry, rx - 1);
      end
    end
  end
  
  % LAMDA IN THE RIGHT REGION
  LAMDAR = zeros(RY, RX);
  for ry = 1: RY
    for rx = RX: -1: 1
      if (rx == RX)
        LAMDAR(ry, rx) = LAMDA0(ry, rx);
      else
        LAMDAR(ry,rx) = LAMDA0(ry, rx + 1);
      end
    end
  end
  
  % LAMDA IN THE BOTTOM REGION
  LAMDAB = zeros(RY, RX);
  for rx = 1: RX
    for ry = 1: RY
      if (ry == 1)
        LAMDAB(ry, rx) = LAMDA0(ry, rx);
      else
        LAMDAB(ry, rx) = LAMDA0(ry - 1, rx);
      end
    end
  end
  
  % LAMDA IN THE TOP REGION
  LAMDAT = zeros(RY, RX);
  for rx = 1: RX
    for ry = RY: -1: 1
      if (ry == RY)
        LAMDAT(ry, rx) = LAMDA0(ry, rx);
      else
        LAMDAT(ry, rx) = LAMDA0(ry + 1, rx);
      end
    end
  end

end