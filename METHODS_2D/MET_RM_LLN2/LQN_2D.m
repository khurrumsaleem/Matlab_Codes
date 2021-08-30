% LEVEL SYMMETRY QUADRATURE FOR 2D IMPLEMENTATIONS
% MAXIMUM QUADRATURE ORDER OF 18 THAT GENERATE 180 DIRECTIONS

% INPUT: N = QUADRATURE ORDER

% OUTPUT: QUAD(miu, ni, w)
%           miu = DISCRETE ORDINATES IN X
%           ni = DISCRETE ORDINATES IN Y
%           w = WEIGHTS
%         chi = DATABASE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QUAD, chi] = LQN_2D(N)

  % VARIABLE INITIALIZATIONS
  M = N*(N+2)/2; % NUMBER OF DIRECTIONS IN THE XY PLANE
  miu = zeros(1, M); ni = zeros(1,M); 
  fi = zeros(1, M); w = zeros(1,M);
  
  % DATABASE
  if N == 2
    chi(1) = 0.5773502692; wlist(1) = 1.0;
  elseif N == 4
    chi(1) = 0.3500212; wlist(1) = 0.33333333;
    chi(2) = 0.8688903; 
  elseif N == 6
    chi(1) = 0.2666355; wlist(1) = 0.1761263;
    chi(2) = 0.6815076; wlist(2) = 0.1572071;
    chi(3) = 0.9261808; 
  elseif N == 8
    chi(1) = 0.2182179; wlist(1) = 0.1209877;
    chi(2) = 0.5773503; wlist(2) = 0.0907407;
    chi(3) = 0.7867958; wlist(3) = 0.0925926;
    chi(4) = 0.9511897; 
  elseif N == 10
    chi(1) = 0.1893213; wlist(1) = 0.0893031;
    chi(2) = 0.5088818; wlist(2) = 0.0725292;
    chi(3) = 0.6943189; wlist(3) = 0.0450438;
    chi(4) = 0.8397600; wlist(4) = 0.0539281;
    chi(5) = 0.9634910; 
  elseif N == 12
    chi(1) = 0.1672126; wlist(1) = 0.0707626;
    chi(2) = 0.4595476; wlist(2) = 0.0558811;
    chi(3) = 0.6280191; wlist(3) = 0.0373377;
    chi(4) = 0.7600210; wlist(4) = 0.0502819;
    chi(5) = 0.8722706; wlist(5) = 0.0258513;
    chi(6) = 0.9716377;
  elseif N == 14
    chi(1) = 0.1519859; wlist(1) = 0.0579970;
    chi(2) = 0.4221570; wlist(2) = 0.0489008;
    chi(3) = 0.5773503; wlist(3) = 0.0221497;
    chi(4) = 0.6988921; wlist(4) = 0.0407009;
    chi(5) = 0.8022263; wlist(5) = 0.0393867;
    chi(6) = 0.8936911; wlist(6) = 0.0245518;
    chi(7) = 0.9766272; wlist(7) = 0.0121325;
  elseif N == 16
    chi(1) = 0.1389568; wlist(1) = 0.0489872;
    chi(2) = 0.3922893; wlist(2) = 0.0413296;
    chi(3) = 0.5370966; wlist(3) = 0.0212326;
    chi(4) = 0.6504264; wlist(4) = 0.0256207;
    chi(5) = 0.7467506; wlist(5) = 0.0360486;
    chi(6) = 0.8319966; wlist(6) = 0.0144589;
    chi(7) = 0.9092855; wlist(7) = 0.0344958;
    chi(8) = 0.9805009; wlist(8) = 0.0085179;
  elseif N == 18
    chi(1) = 0.1293445; wlist(1) = 0.0422646;
    chi(2) = 0.3680438; wlist(2) = 0.0376127;
    chi(3) = 0.5041652; wlist(3) = 0.0066907;
    chi(4) = 0.6106625; wlist(4) = 0.0391919;
    chi(5) = 0.7011669; wlist(5) = 0.0042550;
    chi(6) = 0.7812562; wlist(6) = 0.0423662;
    chi(7) = 0.8538662; wlist(7) = 0.0092396;
    chi(8) = 0.9207680; wlist(8) = 0.0156648;
    chi(9) = 0.9831277; wlist(9) = 0.0136576;
                        wlist(10) = 0.0139903;
  end
  
  % FILLING THE ORDINATES
  d = 0; nlevel = length(chi);
  for n = 1: nlevel
    aux = 0;
    for m = nlevel - n + 1: -1: 1
      d = d + 1; aux = aux + 1;
      miu(d) = chi(m); ni(d) = chi(aux); fi(d) = chi(n);
    end
  end
  
  % FILLING THE WEIGHTS
  p = 0;
  for n = 1: M / 4
    if w(n) == 0
      p = p + 1;
      w(n) = wlist(p);
      aux1 = miu(n); aux2 = ni(n); aux3 = fi(n);
      for m = 1: M / 4
        if aux1 == ni(m) && aux2 == miu(m) && aux3 == fi(m)
          w(m) = wlist(p);
        elseif aux1 == ni(m) && aux2 == fi(m) && aux3 == miu(m)
          w(m) = wlist(p);
        elseif aux1 == fi(m) && aux2 == ni(m) && aux3 == miu(m)
          w(m) = wlist(p);
        elseif aux1 == miu(m) && aux2 == fi(m) && aux3 == ni(m)
          w(m) = wlist(p);
        elseif aux1 == fi(m) && aux2 == miu(m) && aux3 == ni(m)
          w(m) = wlist(p);
        end
      end
    end 
  end
  
  % FILLING THE REMAINDER QUADRANTS
  aux = 0;
  for q = 1:4
    for m = 1: M / 4
      aux = aux + 1;
      if q == 2
        miu(aux) = -miu(m); ni(aux) = ni(m);
        fi(aux) = fi(m); w(aux) = w(m);
      elseif q == 3
        miu(aux) = - miu(m); ni(aux) = - ni(m);
        fi(aux) = fi(m); w(aux) = w(m);
      elseif q == 4
        miu(aux) = miu(m); ni(aux) = - ni(m);
        fi(aux) = fi(m); w(aux) = w(m);
      end
    end
  end
  
  QUAD = [miu', ni', w'];
  
end