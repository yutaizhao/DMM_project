
/*** Input Parameters ***/

L = 1.0;          // Length of the bar
E = 1e4;          // Young'A Modulus
A = 1.0;          // Surface
Fd = 500;         // Applied force at the last element (N)
N = 100;           // Number of elements
n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

eles = 5;           // Number of elements in each subdomain 
nodes = eles + 1;   // Number of nodes in each subdomain
S = N / eles;      // Number of subdomains
Inodes = S + 1;        // Number of interfacial nodes
H = L / S;         // Length of the subdomains




/***********************************
* Section 1 : Resolution without DDM
***********************************/
 
// Stiffness matrix for one element (local matrix)
k_local = (E * A / h) * [1, -1; -1, 1];

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
disp("Nodal Deplacements by elimination:");
disp(u);
//disp("Nodal Forces:");
//disp(f);

//Penality method
M = (E * A / h) * 1e4;          //penality value
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
//disp("Nodal Force by penalty:");
//disp(f_penal);


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
//disp("Nodal Forces by lagrangian:");
//disp(f);



/*********************************
* Section 2 : Domain Decomposition
* Primal Schur Approach
**********************************/

// Initialization of stiffness matices of subdomain 1
K_1 = zeros(nodes, nodes);
K_1_strim = zeros(nodes-1, nodes-1);
for i = 1:eles
    K_1(i:i+1, i:i+1) = K_1(i:i+1, i:i+1) + k_local;
end
K_1_strim = K_1(2:$,2:$);

// Initialization of stiffness matices of subdomain s
K_s = zeros(nodes, nodes);
for i = 1:eles-1
    K_s(i, i) = K_s(i, i) + 2; 
end
for i = 1:eles-2
    K_s(i, i+1) = K_s(i, i+1) - 1; 
    K_s(i+1, i) = K_s(i+1, i) - 1;
end
K_s(1, eles) = K_s(1, eles) - 1;
K_s(eles, 1) = K_s(eles, 1) -1;
K_s(nodes, eles-1) = K_s(nodes, eles-1) - 1;
K_s(eles-1, nodes) = K_s(eles-1, nodes) - 1;
K_s(eles, eles) = K_s(eles, eles) + 1;
K_s(nodes, nodes) = K_s(nodes, nodes) + 1; 
K_s = K_s * (E * A / h);

// Initialization of stiffness matices of subdomain S
K_S = zeros(nodes, nodes);
for i = 1 : eles-1
    K_S(i, i) = K_S(i, i) + 2;
    K_S(i, i+1) = K_S(i, i+1) - 1;
    K_S(i+1, i) = K_S(i+1, i) - 1;
end
K_S(eles, eles) = K_S(eles, eles) + 1;
K_S(nodes, nodes) = K_S(nodes, nodes) + 1;
K_S(1, nodes) = K_S(1, nodes) - 1;
K_S(nodes, 1) = K_S(nodes, 1) - 1;
K_S = K_S * (E * A / h);

//disp("init K_1:");
//disp(K_1);
//disp("K_1_strim:");
//disp(K_1_strim);
//disp("init K_s:");
//disp(K_s);
//disp("init K_S:");
//disp(K_S);

// Get sizes 
[k1s_r,k1s_c] = size(K_1_strim);
[ks_r,ks_c] = size(K_s);
[kS_r,kS_c] = size(K_S);

// Calculate Primal Schur Matrices of subdomains
Sp_1 = K_1_strim(k1s_r, k1s_c)        -K_1_strim(k1s_r, 1:k1s_r-1) * inv(K_1_strim(1:k1s_r-1, 1:k1s_r-1)) * K_1_strim(1:k1s_r-1, k1s_c);
Sp_s = K_s(ks_r-1:ks_r, ks_c-1:ks_c)  -K_s(ks_r-1:ks_r, 1:ks_r-2)  * inv(K_s(1:ks_r-2, 1:ks_r-2))         * K_s(1:ks_r-2,  ks_c-1:ks_c);
Sp_S = K_S(kS_r, kS_c)                -K_S(kS_r       , 1:kS_r-1)  * inv(K_S(1:kS_r-1, 1:kS_r-1))         * K_S(1:kS_r-1,  kS_c); 

//disp("Sp_1:");
//disp(Sp_1);
//disp("Sp_s:");
//disp(Sp_s);
//disp("Sp_S:");
//disp(Sp_S);

// The concatenated Primal complements :

Sp = zeros((2*S-2,2*S-2));
for i = 1:S-2
    Sp(2*i:2*i+1,2*i:2*i+1)= Sp_s;
end
Sp(1,1)=Sp_1;
Sp(2*S-2,2*S-2)=Sp_S;

//disp("Sp:");
//disp(Sp);

// Calculate Primal assembly operators of subdomains

A = zeros((S-1,2*S-2));
for i = 2:S-1
    A(i-1,2*i-2) = 1;
    A(i,2*i-1) = 1;
end
A(1,1)=1;
A(S-1,2*S-2)=1;

//disp("A:");
//disp(A);

// Calculate Primal schur second members of subdomains

bp = zeros((2*S-2,1));
bp(2*S-2) = Fd;

//disp("bp:");
//disp(bp);

// CCL
SSp = A*Sp*A'
BBp = A*bp

Up = SSp\BBp

disp("Subdomains Deplacements by Primal:");
disp(Up);

/*********************************
* Section 2 : Domain Decomposition
* Dual Schur Approach
**********************************/

// Calculate Dual assembly operators of subdomains

A_ = zeros((S-1,2*S-2));
for i = 2:S-1
    A_(i-1,2*i-2) = -1;
    A_(i,2*i-1) = 1;
end
A_(1,1)=1;
A_(S-1,2*S-2)=-1;

//disp("A_:");
//disp(A_);

// The concatenated Dual complements :

Sd_1 = pinv(Sp_1)
Sd_s = pinv(Sp_s)
Sd_S = pinv(Sp_S)

Sd = zeros((2*S-2,2*S-2));
for i = 1:S-2
    Sd(2*i:2*i+1,2*i:2*i+1)= Sd_s;
end
Sd(1,1)=Sd_1;
Sd(2*S-2,2*S-2)=Sd_S;

//disp("Sd:");
//disp(Sd);

//The Concatenated rigid body modes 

Rb_1 = 0; //eliminated
Rb_s = ones(2,1); //kernel
Rb_S = 1;

Rb = zeros((2*S-2,S-1));

for i = 1:S-2
    Rb(2*i,i)= Rb_s(1);
    Rb(2*i+1,i)= Rb_s(2);
end
Rb(2*S-2,S-1) = Rb_S;

//disp("Rb:");
//disp(Rb);

// Dual formulation of the interface problem 

SSd = A_ * Sd * A_';
G = A_ * Rb;
bbd = A_ * Sd * bp ;
e = Rb' * bp;

//disp("SSd:");
//disp(SSd);

//disp("G:");
//disp(G);

//disp("bbd:");
//disp(bbd);

//disp("e:");
//disp(e);

[SSd_r, SSd_c] = size(SSd);
[G_r, G_c] = size(G);
[Gt_r, Gt_c] = size(G');
[bbd_r, bbd_c] = size(bbd);
[e_r, e_c] = size(e);

// CCL

dual_M = zeros((SSd_r+Gt_r, SSd_c+G_c));
dual_rhs = zeros((bbd_r+e_r,1));

dual_M(1:SSd_r, 1:SSd_c)= SSd;
dual_M(1:SSd_r, SSd_c+1:$)= G;
dual_M(SSd_r+1:$, 1:SSd_c)= G';

dual_rhs(1:bbd_r) = -bbd;
dual_rhs(bbd_r+1:$) = -e;

disp("dual_M:");
disp(dual_M);

disp("dual_rhs:");
disp(dual_rhs);

dual_Sol = dual_M\dual_rhs
Ud = dual_Sol(bbd_r+1:$)

disp("Subdomains Deplacements by Dual:");
disp(Ud);

