% HOMOGENEOUS MODEL PROBLEM (Dominguez et. el., 2010)

  % QUADRATURE ORDER
  N = 4;
  
  % ZONES (1. st; 2. ss)
  ZON = [1; 
         0.95];
  
  % REGIONS IN X (1. xlength; 2. xcells)
  XDOM = [10, 10; 
          10, 10];
  
  % REGIONS IN Y (1. ylength; 2. ycells)
  YDOM = [10, 10; 
          10, 10];
  
  % ZONE MAP Y x X ( # zona ) 
  ZMAP = [1, 1;
          1, 1];
  
  % SOURCE MAP Y x X ( Q ) 
  QMAP = [1, 1;
          1, 1];
  
  % BOUNDARY CONDITIONS (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [0.0 0.0 0.0 0.0];
  
  % TOLERANCE
  TOL = 0.00001;