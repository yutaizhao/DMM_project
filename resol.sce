// Input Parameters
L = 1.0;          // Length of the bar
E = 1e4;           // Young's Modulus
S = 1;            // 
N = 10;           // Number of elements
h = L / N;        // distance
n = N + 1;        // Number of nodes
Fd = 500;         // Applied force at the end (N)


// Stiffness matrix for one element (local matrix)
k_local = (E * S / h) * [1, -1; -1, 1];

// Global Stiffness Matrix
K_global = zeros(n, n);

for i = 1:n-1
    K_global(i:i+1, i:i+1) = K_global(i:i+1, i:i+1) + k_local;
end

//Direct elimination of DOF
K_trim = K_global(2:$, 2:$); // Trim first row and column, cuz u[0] = 0
f_trim = zeros(n-1, 1); // Trim first value, cuz u[0] = 0
f_trim($) = Fd ; // Given force
u_trim = K_trim \ f_trim; //Solve
u = [0; u_trim]; //Solution of u
f = K_global*u; //Solution od f

disp("Stiffness Matrix:");
//disp(K_global);
disp("Nodal Deplacements:");
disp(u);
disp("Nodal Forces:");
disp(f);

//Penality method
M = (E * S / h) * 1e4;          //penality value
K_penal = K_global;
K_penal(1,1) = K_penal(1,1) + M;
f_penal = zeros(n, 1);
//f_penal(1) = 1; // set initial guess to f1 
f_penal($) = Fd;
u = K_penal \ f_penal;

f_penal(1) = -M * u(1); // re-calculate f1 (final solution)

disp("Stiffness Matrix:");
disp(K_penal);
disp("Nodal Deplacements:");
disp(u);
disp("Nodal Forces:");
disp(f_penal);


//Langrangian multiplier
K_extend = zeros(n+1, n+1);
K_extend(1:n,1:n) = K_global;
K_extend(1,n+1)=1;
K_extend(n+1,1)=1;
f_lm = zeros(n+1,1);
f_lm($-1) = Fd; // given value for right side

x_extend = K_extend \ f_lm;

u = x_extend(1:$-1);
lambda = x_extend($);
f = f_lm(1:$-1);
f(1) = -lambda;

disp("Stiffness Matrix:");
//disp(K_extend);
disp("Nodal Deplacements:");
disp(u);
disp("Nodal Forces:");
disp(f);
