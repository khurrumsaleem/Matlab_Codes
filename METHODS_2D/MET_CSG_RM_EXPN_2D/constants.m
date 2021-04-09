% CALCULATE SOME AUXILIARY CONSTANTS USED IN THE RESPONSE MATRIX FUNCTION

function [XD0, XDB, XDT, YD0, YDL, YDR, XN0, XNB, XNT, YN0, YNL, YNR] = ...
          constants(QUAD, ZON, XDOM, YDOM, ZMAP, LAMDA0, LAMDAB, LAMDAT, LAMDAL, LAMDAR)
      
 % SHORTHAND
 % XD0(ry, rx) XDB(ry, rx) XDT(ry, rx) YD0(ry, rx) YDL(ry, rx) YDR(ry, rx)
 % XN0(ry, rx) XNB(ry, rx) XNT(ry, rx) YD0(ry, rx) YNL(ry, rx) YNR(ry, rx)
 
 M = length(QUAD(:,1));
 RX = length(XDOM(1, :)); RY = length(YDOM(1,:));
 
 % CALCULATE THE XDs and YDs CONSTANTS
 XD0 = zeros(RY, RX); XDB = zeros(RY, RX); XDT = zeros(RY, RX);
 YD0 = zeros(RY, RX); YDL = zeros(RY, RX); YDR = zeros(RY, RX);
 for ry = 1: RY
   for rx = 1: RX
     z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z);
     aux0X = 0; auxB = 0; auxT = 0; 
     aux0Y = 0; auxL = 0; auxR = 0;
     for k = 1: M
        miu = QUAD(k,1); theta = QUAD(k, 2); w = QUAD(k,3);
            
        aux0X = aux0X + w / (st + miu * LAMDA0(ry, rx));
        auxB = auxB + w / (st + miu * LAMDAB(ry, rx));
        auxT = auxT + w / (st + miu * LAMDAT(ry, rx));
            
        aux0Y = aux0Y + w / (st + theta * LAMDA0(ry, rx));
        auxL = auxL + w / (st + theta * LAMDAL(ry, rx));
        auxR = auxR + w / (st + theta * LAMDAR(ry, rx));
      end
      aux0X = 0.25 * ss * aux0X; auxB = 0.25 * ss * auxB; 
      auxT = 0.25 * ss * auxT;   aux0Y = 0.25 * ss * aux0Y;
      auxL = 0.25 * ss * auxL;   auxR = 0.25 * ss * auxR;
         
      XD0(ry, rx) = 1 - aux0X;  YD0(ry, rx) = 1 - aux0Y;
      XDB(ry, rx) = 1 - auxB;   YDL(ry, rx) = 1 - auxL;
      XDT(ry, rx) = 1 - auxT;   YDR(ry, rx) = 1 - auxR;
   end
 end

 % CALCULATE THE Ns CONSTANTS
 XN0 = zeros(RY, RX); XNB = zeros(RY, RX); XNT = zeros(RY, RX); 
 YN0 = zeros(RY, RX); YNL = zeros(RY, RX); YNR = zeros(RY, RX);
 for ry = 1: RY
   for rx = 1: RX
     H = XDOM(1, rx); K = YDOM(1, ry);
     XN0(ry, rx) = 0.5 * H * LAMDA0(ry, rx) / sinh(0.5 * H * LAMDA0(ry, rx));
     XNB(ry, rx) = 0.5 * H * LAMDAB(ry, rx) / sinh(0.5 * H * LAMDAB(ry, rx));
     XNT(ry, rx) = 0.5 * H * LAMDAT(ry, rx) / sinh(0.5 * H * LAMDAT(ry, rx));
         
     YN0(ry, rx) = 0.5 * K * LAMDA0(ry, rx) / sinh(0.5 * K * LAMDA0(ry, rx));
     YNL(ry, rx) = 0.5 * K * LAMDAL(ry, rx) / sinh(0.5 * K * LAMDAL(ry, rx));
     YNR(ry, rx) = 0.5 * K * LAMDAR(ry, rx) / sinh(0.5 * K * LAMDAR(ry, rx));
   end
 end
 
end