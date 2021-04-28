% CALCULATE SOME AUXILIARY CONSTANTS USED IN THE RESPONSE MATRIX FUNCTION 
% VERSION 4

function [XD0, XDB, XDT, YD0, YDL, YDR, XN0, XNB, XNT, YN0, YNL, YNR] = ...
          constants(QUAD, ZON, XDOM, YDOM, ZMAP, LANDA0, LANDAB, LANDAT, LANDAL, LANDAR)
      
 % SHORTHAND
 % XD0(JB, IB) XDB(JB, IB) XDT(JB, IB) YD0(JB, IB) YDL(JB, IB) YDR(JB, IB)
 % XN0(JB, IB) XNB(JB, IB) XNT(JB, IB) YD0(Jb, IB) YNL(JB, IB) YNR(JB, IB)
 
 M = length(QUAD(:,1));
 RX = length(XDOM(1, :)); RY = length(YDOM(1,:));
 I = sum(XDOM(2, :));     J = sum(YDOM(2, :));
 
 % CALCULATE THE XDs and YDs CONSTANTS
 XD0 = zeros(J, I); XDB = zeros(J, I); XDT = zeros(J, I);
 YD0 = zeros(J, I); YDL = zeros(J, I); YDR = zeros(J, I);
 JB = 0;
 for ry = 1: RY
   ncy = YDOM(2, ry);
   for j = 1: ncy
      JB = JB + 1;
      IB = 0;
      for rx = 1: RX
        z = ZMAP(ry, rx); st = ZON(1, z); ss = ZON(2, z);
        ncx = XDOM(2, rx);
        for i = 1: ncx
          IB = IB + 1;
          aux0X = 0; auxB = 0; auxT = 0; 
          aux0Y = 0; auxL = 0; auxR = 0;
          for k = 1: M
            miu = QUAD(k,1); theta = QUAD(k, 2); w = QUAD(k,3);
            
            aux0X = aux0X + w / (st + miu * LANDA0(JB, IB));
            auxB = auxB + w / (st + miu * LANDAB(JB, IB));
            auxT = auxT + w / (st + miu * LANDAT(JB, IB));
            
            aux0Y = aux0Y + w / (st + theta * LANDA0(JB, IB));
            auxL = auxL + w / (st + theta * LANDAL(JB, IB));
            auxR = auxR + w / (st + theta * LANDAR(JB, IB));
          end
          aux0X = 0.25 * ss * aux0X; auxB = 0.25 * ss * auxB; 
          auxT = 0.25 * ss * auxT;   aux0Y = 0.25 * ss * aux0Y;
          auxL = 0.25 * ss * auxL;   auxR = 0.25 * ss * auxR;
         
          XD0(JB, IB) = 1 - aux0X;  YD0(JB, IB) = 1 - aux0Y;
          XDB(JB, IB) = 1 - auxB;   YDL(JB, IB) = 1 - auxL;
          XDT(JB, IB) = 1 - auxT;   YDR(JB, IB) = 1 - auxR;
       end
     end
   end
 end
 
 % CALCULATE THE Ns CONSTANTS
 XN0 = zeros(J, I); XNB = zeros(J, I); XNT = zeros(J, I); 
 YN0 = zeros(J, I); YNL = zeros(J, I); YNR = zeros(J, I);
 JB = 0;
 for ry = 1: RY
   leny = YDOM(1, ry); ncy = YDOM(2, ry); yh = leny / ncy;
   for j = 1: ncy
     JB = JB + 1;
     IB = 0;
     for rx = 1: RX
       lenx = XDOM(1, rx); ncx = XDOM(2, rx); xh = lenx / ncx;
       for i = 1: ncx
         IB = IB + 1;
         XN0(JB, IB) = 1;
         XNB(JB, IB) = 1;
         XNT(JB, IB) = 1;
         
         YN0(JB, IB) = 1;
         YNL(JB, IB) = 1;
         YNR(JB, IB) = 1;
       end
     end
   end
 end
 
end