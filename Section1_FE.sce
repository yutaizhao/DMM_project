// direct solution without domain decomposition 
function K = K_global(E, A, h, n)
    K = zeros(n, n);
    K = 2 * eye(n, n) + diag(-ones(n - 1, 1), 1) + diag(-ones(n - 1, 1), -1);
    K(1, 1) = 1;
    K(n, n) = 1;
    K = K * E * A / h;
endfunction

// direct elimination
function [u, f] = elimination(n, K_global, Fd)
    K_trim = K_global(2:$, 2:$);
    f_trim = zeros(n - 1, 1);
    f_trim($) = Fd;
    u_trim = K_trim \ f_trim; 
    u = [0; u_trim];
    f = K_global * u;
endfunction

// penalty method
function [u, f] = penalty(E, A, h, K_global, Fd)
    M = (E * A / h) * 1e4;
    K_penal = K_global;
    K_penal(1,1) = K_penal(1,1) + M;
    f = zeros(n, 1);
    f($) = Fd;
    u = K_penal \ f;
    f(1) = - M * u(1);
endfunction    

// lagrangian multiplier
function [u, f] = lagrangian(n, K_global, Fd)
    K_extend = zeros(n + 1, n + 1);
    K_extend(1:n, 1:n) = K_global;
    K_extend(1 ,n + 1) = 1;
    K_extend(n + 1, 1) = 1;
    f_lm = zeros(n + 1, 1);
    f_lm($ - 1) = Fd;
    x_extend = K_extend \ f_lm;
    u = x_extend(1:$ - 1);
    f = f_lm(1:$ - 1);
    f(1) = - x_extend($);
endfunction
