exec("Section1_FE.sce");
exec("direct_solvers.sce");
exec("iterative_solvers.sce");
exec("Tools.sce");

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
eles = 5;

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

// extract u_ref for interfacial nodes
function u_ref = extract_u_ref(u, S, eles)
    u_ref = zeros(S-1, 1);
    for i = 1:S-1
        u_ref(i) = u(i * eles + 1);
    end
endfunction


tol = 1e-6;
max_iter = 1000;
K = K_global(E, A, h, n);

/* Executions */ 

/* Section1 */ 
[u, f] = elimination(n, K, Fd);
[u_pen, f_pen] = penalty(E, A, h, K, Fd);
[u_L, f_L] =lagrangian(n, K, Fd);

u_ref = extract_u_ref(u, S, eles);

/* Section2,3 - Direct */ 
Up=Primal_direct(eles,S,E,A,h,Fd);
Ud=Dual_direct(eles,S,E,A,h,Fd);

/* Section2 - Primal CG */ 
uub=Primal_Conjugate_Gradient(eles,S,E,A,h,Fd,max_iter,tol)

/* Section3 - Primal CG precond */
u_BDD = Primal_BDD(eles, S, E, A, h, Fd, max_iter, tol);

/* Section2 - Dual CG precond */ 
u_b_conca_extract = Dual_TEFI(eles,S,E,A,h,Fd,max_iter,tol); 

disp("u ref:");
disp(u);
disp(u_pen);
disp(u_L);

disp("u Primal:");
disp(Up);
disp("u Dual:");
disp(Ud);

disp("uub:");
disp(uub);
disp("u_BDD:");
disp(u_BDD);
disp("u_b_conca_extract:");
disp(u_b_conca_extract);

