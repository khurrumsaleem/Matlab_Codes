% REF: Mello et. al., 2002

  % QUADRATURE ORDER
  N = 6;
  
  % ZONES (1. st; 2. ss)
  ZON = [0.330263, 0.499122, 0.694676; 
          0.314419, 0.494460, 0.634833];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [8, 16, 8, 24; 
          8, 16, 8, 24];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [8, 8, 8, 24, 8, 8; 
          8, 8, 8, 24, 8, 8];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 3, 2, 1;
          1, 3, 2, 1;
          1, 3, 2, 1;
          1, 3, 2, 1;
          1, 3, 2, 1;
          1, 3, 2, 1];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [0.0, 0.0, 1.0, 0.0;
          0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0;
          0.0, 0.0, 0.0, 0.0];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [0.0 -1.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.00001;