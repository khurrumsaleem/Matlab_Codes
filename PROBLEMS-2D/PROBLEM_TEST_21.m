% SYMMETRIC PROBLEM WITH 4 REGIONS

  % QUADRATURE ORDER
  N = 4;
  
  % ZONES (1. st; 2. ss)
  ZON = [1; 
         0.95];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [1, 1; 
          10, 10];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [1, 1; 
          10, 10];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 1;
          1, 1];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1.0, 1.0;
          1.0, 1.0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [0.0 0.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.00001;