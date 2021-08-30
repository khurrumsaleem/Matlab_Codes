% SHORTHAND RM_EXPN_2D

% METHOD CALL
[SCALAR_FLUX, ...  % SCALAR FLUX IN EACH NODE
 ANGULAR_FLUX, ... % ANGULAR FLUX IN EACH NODE
 X_ANG_FLUX, ...   % ANGULAR FLUX AT X EDGES
 Y_ANG_FLUX, ...   % ANGULAR FLUX AT Y EDGES
 ITER, ...         % NUMBER OF ITERATIONS
 TIME ...          % CPU TIME
 ] = MET_RM_LLN_2D_3(N, ZON, XDOM, YDOM, ZMAP, QMAP, BC, TOL);

% POST-PROCESSING
POST_2D(N, ZON, XDOM, YDOM, ZMAP, BC, SCALAR_FLUX, X_ANG_FLUX, Y_ANG_FLUX);