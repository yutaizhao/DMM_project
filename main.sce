exec("direct_solvers.sci");
exec("functions.sci");

/*** Input Parameters ***/

// fixed values
E = 1e4;                   // Young'A Modulus
A = 1.0;                   // Surface
Fd = 500;                  // Applied force at the last element (N)


// define total number of elements and number of elements in each subdomain
N_list = [60, 120, 180, 240, 300];
eles_list = [2, 3, 4, 5, 10, 15, 20, 25];

// define length of the bar
L_list = [1.0, 2.0, 3.0, 4.0, 5.0];

// test values
L = 1.0;
N = 15;
eles = 3;

n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

nodes = eles + 1;   // Number of nodes in each subdomain
S = N / eles;      // Number of subdomains
Inodes = S + 1;        // Number of interfacial nodes
H = L / S;         // Length of the subdomains

// relative error of two vectors
function err = relative_error(x_ref, x)
    err = norm(x - x_ref) / norm(x_ref);
endfunction

K = K_global(E, A, h, n);
[u, f] = elimination(n, K, Fd);
[u_pen, f_pen] = penalty(E, A, h, K, Fd);
[u_L, f_L] =lagrangian(n, K, Fd);

Up = Conjugate_Gradient_Dual_Schur(eles, S, E, A, h, Fd, 100, 0.0001);

//err = relative_error(u, Up);

disp("global K");
disp(K);
disp("u ref:");
disp(u);
disp(u_pen);
disp(u_L);
disp("u precond schur:");
disp(Up);
//disp("error:");
//disp(err);
