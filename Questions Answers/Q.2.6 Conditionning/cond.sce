exec("./Tools.sce");

// fixed values
E = 1e4;                   // Young'A Modulus
A = 1.0;                   // Surface
Fd = 500;                  // Applied force at the last element (N)
L=20;

N = 90;
n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

eles = 30;


S = N / eles;      // Number of subdomains
nodes = eles + 1;   // Number of nodes in each subdomain
Inodes = S + 1;        // Number of interfacial nodes
H = L / S;         // Length of the subdomains

SSp_inv=build_SSp_inv(S);
Sp_conca=build_Sp_conca(S);
AA=build_AA(S);
SSp = AA*Sp_conca*AA';
disp("eles = " + string(eles));
disp(cond(SSp_inv*SSp));


    

