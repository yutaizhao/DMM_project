// Input Parameters
L = 1.0;          // Length of the bar
E = 1e4;           // Young's Modulus
S = 1.0;            // 
N = 100;           // Number of elements
h = L / N;        // distance
n = N + 1;        // Number of nodes
Fd = 500;         // Applied force at the end (N)

eles = 5;         // Number of elements in each subdomain 
nodes = eles + 1;    // Number of nodes in each subdomain
ND = N / eles;;        // Number of subdomains
H = L / ND;
nd = ND + 1;        // Number of interfacial nodes

 
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

//disp("Stiffness Matrix:");
//disp(K_global);
//disp("Nodal Deplacements:");
//disp(u);
//disp("Nodal Forces:");
//disp(f);

//Penality method
M = (E * S / h) * 1e4;          //penality value
K_penal = K_global;
K_penal(1,1) = K_penal(1,1) + M;
f_penal = zeros(n, 1);
//f_penal(1) = 1; // set initial guess to f1 
f_penal($) = Fd;
u = K_penal \ f_penal;

f_penal(1) = -M * u(1); // re-calculate f1 (final solution)

//disp("Stiffness Matrix:");
//disp(K_penal);
disp("Nodal Deplacements by penalty:");
disp(u);
disp("Nodal Force by penalty:");
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

//disp("Stiffness Matrix:");
//disp(K_extend);
disp("Nodal Deplacements by lagrangian:");
disp(u);
disp("Nodal Forces by lagrangian:");
disp(f);


// Domain Decomposition
// Primal Schur Approach

// Initialization of stiffness matices of subdomains
K_1 = zeros(nodes, nodes);
for i = 1:eles
    K_1(i:i+1, i:i+1) = K_1(i:i+1, i:i+1) + k_local;
end


K_d = zeros(nodes, nodes);
for i = 1:eles-1
    K_d(i, i) = K_d(i, i) + 2; 
end
for i = 1:eles-2
    K_d(i, i+1) = K_d(i, i+1) - 1; 
    K_d(i+1, i) = K_d(i+1, i) - 1;
end
K_d(1, eles) = K_d(1, eles) - 1;
K_d(eles, 1) = K_d(eles, 1) -1;
K_d(nodes, eles-1) = K_d(nodes, eles-1) - 1;
K_d(eles-1, nodes) = K_d(eles-1, nodes) - 1;
K_d(eles, eles) = K_d(eles, eles) + 1;
K_d(nodes, nodes) = K_d(nodes, nodes) + 1; 
K_d = K_d * (E * S / h);


K_ND = zeros(nodes, nodes);
for i = 1 : eles-1
    K_ND(i, i) = K_ND(i, i) + 2;
    K_ND(i, i+1) = K_ND(i, i+1) - 1;
    K_ND(i+1, i) = K_ND(i+1, i) - 1;
end
K_ND(eles, eles) = K_ND(eles, eles) + 1;
K_ND(nodes, nodes) = K_ND(nodes, nodes) + 1;
K_ND(1, nodes) = K_ND(1, nodes) - 1;
K_ND(nodes, 1) = K_ND(nodes, 1) - 1;
K_ND = K_ND * (E * S / h);

disp("init K_1:");
disp(K_1);
disp("init K_d:");
disp(K_d);
disp("init K_ND:")
disp(K_ND);


// Calculate Primal Schur Matrices of subdomains
Sp_1 = K_1(nodes, nodes) - K_1(nodes, 1:eles) * K_1(1:eles, 1:eles) * K_1(1:eles, 1);
Sp_d = K_d(eles:nodes, eles:nodes) - K_d(eles:nodes, 1:eles-1) * K_d(1:eles-1, 1:eles-1) * K_d(1:eles-1, eles:nodes);
Sp_ND = K_ND(nodes, nodes) - K_ND(nodes, 1:eles) * K_ND(1:eles, 1:eles) * K_ND(1:eles, 1); 

disp("Sp_1:");
disp(Sp_1);
disp("Sp_d:");
disp(Sp_d);
disp("Sp_ND:");
disp(Sp_ND);

// Calculate Primal schur second members of subdomains



