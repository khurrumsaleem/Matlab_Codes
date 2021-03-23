% CALCULATE THE GAUSS-LEGENDRE QUADRATURE PARAMETERS

% INPUT: N, QUADRATURE ORDER

% OUTPUT: QUAD = [miu; w]
%                miu = DISCRETE-ORDINATES COORDINATES
%                w = WEIGHTS

function QUAD = GAUSS_QUAD(N)

    % INITIAL SETS
    a = -1; b = 1; N = N - 1; N1 = N + 1; N2 = N + 2;
    xu = linspace(-1, 1, N1)';

    % INITIAL GUESS
    y = cos((2 * (0: N)' + 1) * pi / (2 * N + 2)) + (0.27 / N1) * sin(pi * xu *N / N2);

    % LEGENDRE-GAUSS VANDERMONDE MATRIX (LGVM)
    L = zeros(N1,N2);

    % DERIVATIVE OF LGVM
    Lp = zeros(N1,N2);

    % COMPUTE THE ZEROS OF THE N + 1 LEGENDRE POLYNOMIAL USING THE RECURSION RELATION AND THE NEWTON-RAPHSON METHOD
    y0 = 2;
    while max(abs(y - y0)) > eps
        L(:, 1) = 1;   Lp(:, 1) = 0; L(:, 2) = y;   Lp(:, 2) = 1;
        for k = 2: N1
            L(:, k + 1) = ( (2 * k -  1) * y .* L(:, k) - (k - 1) * L(:, k - 1) ) / k;
        end
        Lp = (N2) * ( L(:, N1) - y .* L(:, N2) ) ./ (1 - y .^ 2);   
        y0 = y;
        y = y0 - L(:, N2) ./ Lp;
    end

    % LINEAR MAP FROM [-1, 1] to [a, b]
    miu = (a * (1 - y) + b * (1 + y)) / 2;      

    % COMPUTE THE WEIGHTS
    w = (b - a) ./ ((1 - y .^ 2) .* Lp .^2 ) * (N2 / N1) ^2;
    
    % CALCULATE NEGATIVE ORDINATES
    N = N + 1;
    for m = 1 : N/2
      miu(N / 2 + m) = - miu(m); w(N / 2 + m) = w(m);
    end

    % FORMAT OUTPUT
    QUAD = [miu, w];
    
end