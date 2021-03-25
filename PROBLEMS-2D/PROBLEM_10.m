% REF: Dominguez et. al., 2010

  % QUADRATURE ORDER
  N = 4;
  
  % ZONES (1. st; 2. ss)
  ZON = [1.0, 2.0; 
         0.5, 0.1];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [5, 5; 
          20, 20];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [5, 5; 
          20, 20];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 2;
          2, 2];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1.0, 0.0;
          0.0, 0.0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [-1.0 -1.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.000001;