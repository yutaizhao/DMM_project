
/*** Input Parameters ***/

L = 1.0;          // Length of the bar
E = 1e4;          // Young'A Modulus
A = 1.0;          // Surface
Fd = 500;         // Applied force at the last element (N)
N = 15;           // Number of elements
n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

eles = 5;           // Number of elements in each subdomain 
nodes = eles + 1;   // Number of nodes in each subdomain
S = N / eles;      // Number of subdomains
Inodes = S + 1;        // Number of interfacial nodes
H = L / S;         // Length of the subdomains

/*********************************
* Section 2 : Domain Decomposition
* Primal Schur Approach - TEFI
**********************************/

/* Local stiffness matrix  */

function K=build_K(eles,E,A,h,s,S)
    
    nodes = eles+1;
    k_local = (E * A / h) * [1, -1; -1, 1];
    
    if s==1 then
        // Initialization of stiffness matices of subdomain 1
        K_1 = zeros(nodes, nodes);
        K = zeros(nodes-1, nodes-1);
        for i = 1:eles
            K_1(i:i+1, i:i+1) = K_1(i:i+1, i:i+1) + k_local;
        end
        K = K_1(2:$,2:$);
        
    elseif s==S then 
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
        K = K_S * (E * A / h);
        
    elseif (s > 1 && s < S) then 
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
        K = K_s * (E * A / h);
    end

endfunction

//Extract different parts of K
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

/* Local Sp */

function Sp=build_Sp(s,S)
    if s==1 then
        K_1 = build_K(eles,E,A,h,s,S);
        Sp = K_1($, $)    -K_1($, 1:$-1) * inv(K_1(1:$-1, 1:$-1)) * K_1(1:$-1, $);
    elseif s==S then 
        K_S = build_K(eles,E,A,h,s,S);
        Sp = K_S($, $)    -K_S($, 1:$-1) * inv(K_S(1:$-1, 1:$-1)) * K_S(1:$-1, $); 
    elseif (s > 1 && s < S) then 
        K_s = build_K(eles,E,A,h,s,S);
        Sp = K_s($-1:$, $-1:$)  -K_s($-1:$, 1:$-2) * inv(K_s(1:$-2, 1:$-2)) * K_s(1:$-2, $-1:$);
    end
endfunction

/* Concatenated Sp */

function Sp_conca=build_Sp_conca(S)
    Sp_conca = zeros((2*S-2,2*S-2));
    Sp_1=build_Sp(1,S);
    Sp_S=build_Sp(S,S);
    
    if S>2 then
        Sp_s=build_Sp(2,S);
        for i = 1:S-2
            Sp_conca(2*i:2*i+1,2*i:2*i+1)= Sp_s;
        end
    end
    Sp_conca(1,1)=Sp_1;
    Sp_conca(2*S-2,2*S-2)=Sp_S;
endfunction


/* Local Sd */
function Sd=build_Sd(s,S)
    Sp=build_Sp(s,S);
    Sd = pinv(Sp);
endfunction

/* Global AA_ */
function AA_=build_AA_(S)
    AA_ = zeros((S-1,2*S-2));
    for i = 2:S-1
        AA_(i-1,2*i-2) = -1;
        AA_(i,2*i-1) = 1;
    end
    AA_(1,1)=1;
    AA_(S-1,2*S-2)=-1;
endfunction

/* Local A_ */
function A_=build_A_(s,S)
    if s==1 then
        A_=1;
    elseif s==S then 
        A_=-1;
    elseif (s > 1 && s < S) then 
        A_ = zeros((S-1,2));
        A_(s-1,1) = -1;
        A_(s,2) = 1;
    end
endfunction

/* Local Rb */
function Rb=build_Rb(s,S)
    if s==1 then
        Rb = 0; //useless
    elseif s==S then 
        Rb = 1;
    elseif (s > 1 && s < S) then 
        Sp_s = build_Sp(s,S);
        Rb = kernel(Sp_s);
    end
endfunction

/*Local bp*/
function bp=build_bp(s,S)
    if s==1 then
        bp = 0;
    elseif s==S then 
        bp = Fd;
    elseif (s > 1 && s < S) then 
        bp = zeros(2,1);
    end
endfunction

/*Local bd*/
function bd=build_bd(s,S)
    bp=build_bp(s,S);
    Sd=build_Sd(s,S);
    bd = Sd*bp;
endfunction


/*Local e*/
function e=build_e(s,S)
    Rb=build_Rb(s,S);
    bp=build_bp(s,S);
    e = Rb'*bp;
endfunction



/*Global G*/
function G=build_G(S)
    
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

    Sp_s = K_s($-1:$, $-1:$)  -K_s($-1:$, 1:$-2) * inv(K_s(1:$-2, 1:$-2)) * K_s(1:$-2, $-1:$);

    //The Concatenated rigid body modes 

    Rb_1 = 0; //eliminated
    Rb_s = kernel(Sp_s); //kernel
    Rb_S = 1;

    Rb = zeros((2*S-2,S-1));

    for i = 1:S-2
        Rb(2*i,i)= Rb_s(1);
        Rb(2*i+1,i)= Rb_s(2);
    end
    Rb(2*S-2,S-1) = Rb_S;
    AA_=build_AA_(S);
    G = AA_ * Rb;
endfunction

G = build_G(S);

/*Global e*/
function ee=build_ee(S,Fd)
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

    Rb_1 = 0; //eliminated
    Sp_s =  K_s($-1:$, $-1:$)  -K_s($-1:$, 1:$-2) * inv(K_s(1:$-2, 1:$-2)) * K_s(1:$-2, $-1:$);
    Rb_s = kernel(Sp_s); //kernel
    Rb_S = 1;

    Rb = zeros((2*S-2,S-1));

    for i = 1:S-2
        Rb(2*i,i)= Rb_s(1);
        Rb(2*i+1,i)= Rb_s(2);
    end
    Rb(2*S-2,S-1) = Rb_S;
    bp = zeros((2*S-2,1));
    bp(2*S-2) = Fd;
    
    ee = Rb'*bp;
endfunction

/*A_tild*/
function AA_tild_=build_AA_tild_(S)
    AA_=build_AA_(S);
    AA_tild_ = pinv(AA_*AA_')*AA_;
endfunction

/*Sd_inv*/
function SSd_inv=build_SSd_inv(S)
    AA_tild_=build_AA_tild_(S);
    Sp_conca=build_Sp_conca(S);
    SSd_inv=AA_tild_*Sp_conca*AA_tild_';
endfunction


/*P*/
P = eye(S-1,S-1) - G*inv(G'*G)*G';

/*lambda0*/
ee=build_ee(S,Fd);
lambda0 =  - G*inv(G'*G)*ee;

//TEFI
function uub=TEFI(eles,S,E,A,h,Fd,m,tol)
    
    /*Initialisation */
    
     AA_=build_AA_(S);
     
    //lambda(0)
    lambda = lambda0;
    //Compute local concatenated lambda
    lambda_concatanated = AA_' * lambda;
    
    /** Compute r0 **/
    
    r_concatanated = [];

    //s=1
    bd=build_bd(1,S);
    Sd=build_Sd(1,S);
    lambda = lambda_concatanated(1);
    r_concatanated = [r_concatanated; -bd-Sd*lambda];
    
    //s>1, s<S
    for s=2:S-1
        bd=build_bd(s,S);
        Sd=build_Sd(s,S);
        lambda = lambda_concatanated(2*(s-1):2*(s-1)+1);
        r_concatanated = [r_concatanated; -bd-Sd*lambda];
    end
    
    //s=S
    bd=build_bd(S,S);
    Sd=build_Sd(S,S);
    lambda = lambda_concatanated($);
    r_concatanated = [r_concatanated; -bd-Sd*lambda];

    rr = P'*(AA_*r_concatanated);
    
    /** Compute z0 **/
    
    SSd_inv=build_SSd_inv(S);
    zz = P*SSd_inv*rr;
    
    disp("rr");
    disp(rr);
    disp("zz");
    disp(zz);
    
endfunction

uub = TEFI(eles,S,E,A,h,Fd,100,0.0001);


