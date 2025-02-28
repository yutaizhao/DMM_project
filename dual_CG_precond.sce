
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
        A_ = zeros((S-1,1));
        A_(1,1)=1;
    elseif s==S then 
        A_ = zeros((S-1,2));
        A_(s-1,2*s-2) = -1;
        A_(s,2*s-1) = 1;
    elseif (s > 1 && s < S) then 
        A_ = zeros((S-1,1));
        A_(S-1,1)=-1;
    end
endfunction

/* Local Rb */
function Rb=build_Rb(s,S)
    if s==1 then
        Rb = 0;
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

/*Local G*/
function G=build_G(s,S)
    A_=build_A_(s,S);
    Rb=build_Rb(s,S);
    G = A_*Rb;
endfunction

/*Local e*/
function e=build_e(s,S)
    Rb=build_Rb(s,S);
    bp=build_bp(s,S);
    e = Rb'*bp;
endfunction

/*Local Q*/
function Q=build_Q(S)
    Q = eye(S-1);// for homogeneuos struct
endfunction

/*Local P*/
function P=build_P(s,S)
    Q = build_Q(S);
    G = build_G(s,S);
    P = eye(S-1) - Q*G*inv((G'*Q*G))*G';
endfunction


/*Local lambda*/
function lambda=build_lambda(s,S)
    Q = build_Q(S);
    G = build_G(s,S);
    e=build_e(s,S);
    lambda = - Q*G*inv((G'*Q*G))*G'*e;
endfunction


//TEFI
function uub=TEFI(eles,S,E,A,h,Fd,m,tol)
    
    /*Initialisation */
    r = [];
    for s=1:S
        disp(s);
        bp=build_bp(s,S);
        P = build_P(s,S);
        lambda = build_lambda(s,S);
        Sp=build_Sp(s,S);
        bd = pinv(Sp) * bp;
        r = [r;P'*(-bd-pinv(Sp)*lambda)];
    end
    r_concatanated = AA_ * r;
    
endfunction

uub = TEFI(eles,S,E,A,h,Fd,100,0.0001);


