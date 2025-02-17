
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


function [K_1_strim, K_s, K_S]=build_K(eles,E,A,h)
    
    nodes = eles+1;
    k_local = (E * A / h) * [1, -1; -1, 1];
    
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

endfunction

function [Sp, k1_r, k1_c, ks_r, ks_c, kS_r, kS_c]=build_Sp(S,K_1,K_s,K_S)
    
    [k1_r,k1_c] = size(K_1);
    [ks_r,ks_c] = size(K_s);
    [kS_r,kS_c] = size(K_S);
    
    Sp_1 = K_1(k1_r, k1_c)                -K_1(k1_r, 1:k1_r-1)         * inv(K_1(1:k1_r-1, 1:k1_r-1))         * K_1(1:k1_r-1, k1_c);
    Sp_s = K_s(ks_r-1:ks_r, ks_c-1:ks_c)  -K_s(ks_r-1:ks_r, 1:ks_r-2)  * inv(K_s(1:ks_r-2, 1:ks_r-2))         * K_s(1:ks_r-2,  ks_c-1:ks_c);
    Sp_S = K_S(kS_r, kS_c)                -K_S(kS_r       , 1:kS_r-1)  * inv(K_S(1:kS_r-1, 1:kS_r-1))         * K_S(1:kS_r-1,  kS_c); 
    
    Sp = zeros((2*S-2,2*S-2));
    for i = 1:S-2
        Sp(2*i:2*i+1,2*i:2*i+1)= Sp_s;
    end
    Sp(1,1)=Sp_1;
    Sp(2*S-2,2*S-2)=Sp_S;
endfunction


function AA=build_AA(S)
    AA = zeros((S-1,2*S-2));
    for i = 2:S-1
        AA(i-1,2*i-2) = 1;
        AA(i,2*i-1) = 1;
    end
    AA(1,1)=1;
    AA(S-1,2*S-2)=1;
endfunction

function bp=build_bp(S,Fd)
    bp = zeros((2*S-2,1));
    bp(2*S-2) = Fd;
endfunction

function [Kii, Kib, Kbi, Kbb] = extract_K(K, s, S)
    if (s == 1 || s == S) then
        Kii = K(1:$-1, 1:$-1);
        Kib = K(1:$-1, $);
        Kbi = K($, 1:$-1);
        Kbb = K($, $);
    elseif (s > 1 && s < S) then 
        Kii = K(1:$-2, 1:$-2);
        Kib = K(1:$-2, $-1:$);
        Kbi = K($-1:$, 1:$-2);
        Kbb = K($-1:$, $-1:$);
    end
endfunction


function Up=Conjugate_Gradient(eles,S,E,A,h,Fd,m)
    
    //preliminary 
    [K_1_strim K_s K_S]=build_K(eles,E,A,h);
    [Sp k1_r k1_c ks_r ks_c kS_r kS_c]=build_Sp(S,K_1_strim,K_s,K_S);
    AA = build_AA(S);
    bp = build_bp(S,Fd);
    SSp = AA*Sp*AA';
    BBp = AA*bp;
    
    //ub0
    ub0 = zeros(S-1,1);
    //Compure concatenated residual rb = zeros(2*S-2,1);
    rb_concatanated = [];
    //Compute local concatenated vector
    ub_concatanated = AA' * ub0;
    
    //s=1
    [Kii, Kib, Kbi, Kbb] = extract_K(K_1_strim, 1, S);
    fi_1=zeros(k1_r-1,1);
    fb_1=0;
    ub0_1 = ub_concatanated(1,1);
    ui0_1 = inv(Kii)*(fi_1-Kib*ub0_1);
    res_1 = K_1_strim*[ui0_1; ub0_1]-[fi_1; fb_1] ;
    rb_concatanated = [rb_concatanated; -res_1($)];
    
    //s=2...S-1
    fi_s=zeros(ks_r-2,1);
    fb_s=zeros(2,1);
    for s=2:S-1
        [Kii, Kib, Kbi, Kbb] = extract_K(K_s, s, S);
        ub0_s = ub_concatanated(s:s+1,1);
        ui0_s = inv(Kii)*(fi_s-Kib*ub0_s);
        res_s = K_s*[ui0_s; ub0_s]-[fi_s; fb_s] ;
        rb_concatanated = [rb_concatanated; -res_s($-1:$)];
    end
     
    //s=S
    [Kii, Kib, Kbi, Kbb] = extract_K(K_S, S, S);
    fi_S=zeros(kS_r-1,1);
    fb_S=Fd;
    ub0_S = ub_concatanated($,1);
    ui0_S = inv(Kii)*(fi_S-Kib*ub0_S);
    res_S = K_S*[ui0_S; ub0_S]-[fi_S; fb_S] ;
    rb_concatanated = [rb_concatanated; -res_S($)];
    
    //Compute global residual
    rrb = AA*rb_concatanated;
    ddbk=rrb0;
    
    for k=0:m
        // Compute local vector
        db_concatanated = AA' * ddbk ;
        
        //s=1
        [Kii, Kib, Kbi, Kbb] = extract_K(K_1_strim, 1, S);
        db_1 = db_concatanated(1,1);
        di_1 = inv(Kii)*Kib*db_1;
        res_1 = K_1_strim*[di_1; db_1] ;
        db_concatanated = [db_concatanated; res_1($)];
                           
        //s=2...S-1
        for s=2:S-1
            [Kii, Kib, Kbi, Kbb] = extract_K(K_s, s, S);
            db_s = ub_concatanated(s:s+1,1);
            di_s = inv(Kii)*Kib*db_s;
            res_s = K_s*[di_s; db_s];
            db_concatanated = [db_concatanated; res_s($-1:$)];
        end
     
        //s=S
        [Kii, Kib, Kbi, Kbb] = extract_K(K_S, S, S);
        db_S = db_concatanated(1,1);
        di_S = inv(Kii)*Kib*db_S;
        res_S = K_S*[di_S; db_S] ;
        db_concatanated = [db_concatanated; res_S($)];
        
        
        //Compute global matrix vector product 
       
    end

endfunction



Up=Conjugate_Gradient(eles,S,E,A,h,Fd);




