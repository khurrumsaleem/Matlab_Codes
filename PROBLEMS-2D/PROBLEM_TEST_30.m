% INFINITE MEDIA IN Y

  % Orden de la cuadratura
  N = 4;
  
  % Zonas (1. st; 2. ss)
  ZON = [1; 
         0.95];
  
  % Regiones X (1. xlength; 2. xcells)
  XDOM = [1; 
          2];
  
  % Regiones Y (1. ylength; 2. ycells)
  YDOM = [10; 
          20];
  
  % Mapeo de Zonas Y x X ( # zona ) 
  ZMAP = 1;
  
  % Mapeo de Zonas Y x X ( Q ) 
  QMAP = 1.0;
  
  % Condiciones de contorno (1. left, 2. dowm, 3. right,
  %                          4. up)
  BC = [-1 0 -1 0];
  
  % Tolerancia
  TOL = 0.0001;