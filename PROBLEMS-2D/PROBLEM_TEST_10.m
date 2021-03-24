% PURE ABSORBER PROBLEM WITH INTERNAL SOURCE

  % QUADRATURE ORDER
  N = 2;
  
  % ZONES (1. st; 2. ss)
  ZON = [1; 
         0];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [1; 
          10];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [1; 
          10];
  
  % ZON MAP Y x X ( # zone number ) 
  ZMAP = 1;
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = 1.0;
  
  % BOUNDARY COMDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [0.0 0.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.00001;