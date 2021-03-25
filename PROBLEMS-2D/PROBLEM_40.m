% HIGH DIFFUSIVE PROBLEM

  % QUADRATURE ORDER
  N = 6;
  
  % ZONES (1. st; 2. ss)
  ZON = [0.8, 1; 
         0.0, 0.95];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [1, 44, 5; 
          4, 176, 20];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [1, 44, 5; 
          4, 176, 20];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 2, 2;
          2, 2, 2;
          2, 2, 2];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1.0, 0.0, 0.0;
          0.0, 0.0, 0.0;
          0.0, 0.0, 0.0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [-1.0 -1.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.00001;