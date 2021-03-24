% REF: Dominguez et. al., 2010

  % QUADRATURE ORDER
  N = 4;
  
  % ZONES (1. st; 2. ss)
  ZON = [1.00, 0.10, 0.30; 
          0.50, 0.01, 0.10];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [40, 10, 20, 20, 10; 
          4, 1, 2, 2, 1];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [40, 10, 20, 20, 10; 
          4, 1, 2, 2, 1];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 2, 3, 2, 2;
          2, 2, 3, 2, 2;
          3, 3, 3, 2, 2;
          2, 2, 2, 2, 2;
          2, 2, 2, 2, 2];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1.0, 0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0, 0.0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [-1.0 -1.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.000001;