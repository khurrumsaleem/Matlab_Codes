% BLINDAGEM TEST 1

  % QUADRATURE ORDER
  N = 2;
  
  % ZONES (1. st; 2. ss)
  ZON = [1.00, 1.0; 
         0.95, 0.5];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [5, 5; 
          5, 5];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [5, 5; 
          5, 5];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 2; 
          2, 2];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1, 0; 
          0, 0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [-1 -1 0 0];
  
  % TOLERANCE
  TOL = 0.000001;