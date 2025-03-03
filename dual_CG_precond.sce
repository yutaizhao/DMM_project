
/*** Input Parameters ***/

L = 1.0;          // Length of the bar
E = 1e4;          // Young'A Modulus
A = 1.0;          // Surface
Fd = 500;         // Applied force at the last element (N)
N = 9;           // Number of elements
n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

eles = 3;           // Number of elements in each subdomain 
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

/* Concatenated Sd */
function Sd_conca=build_Sd_conca(S)
    Sd_conca = zeros((2*S-2,2*S-2));
    Sp_1=build_Sp(1,S);
    Sp_S=build_Sp(S,S);
    Sd_1 = pinv(Sp_1);
    Sd_S = pinv(Sp_S);
    
    if S>2 then
        Sp_s=build_Sp(2,S);
        Sd_s = pinv(Sp_s);
        for i = 1:S-2
            Sd_conca(2*i:2*i+1,2*i:2*i+1)= Sd_s;
        end
    end
    
    Sd_conca(1,1)=Sd_1;
    Sd_conca(2*S-2,2*S-2)=Sd_S;
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
        Rb = ones(2,1);
    end
endfunction


/* Global Rb */
function Rb_conca=build_Rb_conca(S)

    Rb_s = ones(2,1); //kernel
    Rb_S = 1;

    Rb_conca = zeros((2*S-2,S-1));

    for i = 1:S-2
        Rb_conca(2*i,i)= Rb_s(1);
        Rb_conca(2*i+1,i)= Rb_s(2);
    end
    
    Rb_conca(2*S-2,S-1) = Rb_S;
    
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

/*Global bp*/
function bp_conca=build_bp_conca(S,Fd)
    bp_conca = zeros((2*S-2,1));
    bp_conca(2*S-2) = Fd;
endfunction

/*Local bd*/
function bd=build_bd(s,S)
    bp=build_bp(s,S);
    Sd=build_Sd(s,S);
    bd = Sd*bp;
endfunction

/*global bd*/
function bbd=build_bbd(S)
    AA_=build_AA_(S);
    Sd_conca=build_Sd_conca(S);
    bp_conca=build_bp_conca(S,Fd);
    bbd = AA_ * Sd_conca * bp_conca ;
endfunction


/*Local e*/
function e=build_e(s,S)
    Rb=build_Rb(s,S);
    bp=build_bp(s,S);
    e = Rb'*bp;
endfunction



/*Global G*/
function G=build_G(S)
    Rb_conca=build_Rb_conca(S);
    AA_=build_AA_(S);
    G = AA_ * Rb_conca;
endfunction

G = build_G(S);

/*Global e*/
function ee=build_ee(S,Fd)
    
    Rb_conca=build_Rb_conca(S);
    bp = zeros((2*S-2,1));
    bp(2*S-2) = Fd;
    
    ee = Rb_conca'*bp;
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
function u_b=TEFI(eles,S,E,A,h,Fd,m,tol)
    
    /*Initialisation */
    
    AA_=build_AA_(S);
    SSd_inv=build_SSd_inv(S);
    
    //lambda(0)
    Lambda = lambda0;

    //Compute local concatenated lambda
    lambda_concatanated = AA_' * Lambda;
    
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
    
    //Check onvergence 
    norm_r0 = norm(rr);  

    /** Compute z0,d0 **/
    
    zz = P*SSd_inv*rr;
    dd = zz;
    
    p_list=[];
    dd_list=[dd];
    
    /* For loop */
    
    for i=1:m 
        
        // Compute local p
        d_concatanated = AA_' * dd;
        local_Sdxd_concatanated = [];
        
        Sd = build_Sd(1,S);
        d = d_concatanated(1);
        local_Sdxd_concatanated = [local_Sdxd_concatanated, Sd*d];
        for s=2:S-1
            Sd = build_Sd(s,S);
            d = d_concatanated(2*(s-1):2*(s-1)+1);
            local_Sdxd_concatanated = [local_Sdxd_concatanated; Sd*d];
        end 
        Sd = build_Sd(S,S);
        d = d_concatanated($);
        local_Sdxd_concatanated = [local_Sdxd_concatanated; Sd*d];
        //END of computation p
        
        p = P'*AA_*local_Sdxd_concatanated;
        p_list=[p_list;p];
        
        // Update CG iter
        a = (rr'*dd)/(d'*p);
        Lambda = Lambda + a*dd;
       
        rr = rr -a*p;
        zz = P*SSd_inv*rr;
        
        beta_list =[];
        for j = 1:i
            beta_list = [beta_list; -(zz'*p_list(j))/(dd_list(j)'*p_list(j))];
        end 
        
        extra_sum = 0;
        for j = 1:i
            extra_sum = extra_sum + beta_list(j)*dd;
        end
        dd = zz + extra_sum;
        dd_list=[dd_list; dd];
        
         //Check onvergence 
        rel_error = norm(rr) / (norm_r0+E-20); 
        if rel_error < tol then
           disp("Converged at iteration " + string(i));
           break;
        end
        
    end
   
    /* Post processing */
    Sd_conca=build_Sd_conca(S);
    SSd = AA_*Sd_conca*AA_';
    bbd=build_bbd(S);
    bp_conca=build_bp_conca(S,Fd);
    Rb_conca=build_Rb_conca(S);
    alpha_b = inv(G'*G)*G'*(-bbd-SSd*Lambda);
    u_b = Sd_conca*(bp_conca+AA_'*Lambda)+Rb_conca*alpha_b;
    disp(u_b);
    
endfunction

u_b = TEFI(eles,S,E,A,h,Fd,10,0.0001);


