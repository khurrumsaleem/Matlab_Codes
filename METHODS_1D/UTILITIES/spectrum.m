% CALCULATE THE EIGENVALUES AND EIGENVECTORS FOR EACH ZONE

% INPUT: QUAD = QUADRATURE
%        ZON = ZONE PARAMETERS

% OUTPUT: vals(m, z) = EIGENVALUE m FOR ZONE z
%         vects(n, m, z) = EIGENVECTOR COMPONENT n ASSOCIATE TO THE
%                          EIGENVALUE m FOR ZONE z

function [vals, vects] = spectrum(QUAD, ZON)
  
  % VARIABLES
  NZ = length(ZON(1,:)); N = length(QUAD(:,1));
  vals = zeros(N, NZ); vects = zeros(N, N, NZ);
  
  % LOOP ALONG ZON
  for z = 1: NZ
    
    % MATERIAL PARAMETER OF THE CURRENT ZONE
    st = ZON(1, z);    ss = ZON(2, z);    c0 = ss / st;
    
    % EIGENVALUES CALCULATION
    if (c0 ~= 0)
       
      for i = 1: N / 2
        miu_i = QUAD(i, 1);
        if i == N / 2
          miu_f = 10;
        else
          miu_f = QUAD(i + 1, 1);
        end
        h = (miu_f - miu_i) / ( 10 ^ N);
        a = miu_i + h; b = miu_f - h;
        r = myRootFunc(a, b, QUAD, c0, 1e-10);
        vals(i,z) = r;
      end
      
    else
        
      for i = 1: N / 2
        miu_i = QUAD(i, 1);
        vals(i, z) = - miu_i;
      end
      
    end
    
    % ORDERING THE EIGENVALUES
    if (c0 ~= 0), vals(:, z) = sort(vals(:, z), 'descend'); end
    for i = 1: N / 2
      vals(N / 2 + i, z) = - vals(i, z); 
    end
    
    % EIGENVECTORS CALCULATION
    if (c0 ~= 0)
        
      for i = 1: N
        for m = 1: N
          miu = QUAD(m, 1);
          vects(m, i, z) = 0.5 * c0 * vals(i, z) / (miu + vals(i, z));
        end
      end
      
    else
        
      for i = 1: N
        for m = 1: N
          if (m == i)
            vects(m, i, z) = 1.0;
          else
            vects(m, i, z) = 0.0;
          end
        end
      end
      
    end
    
  end

end