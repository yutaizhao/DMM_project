
//Create stiffness matrix
function K_s=build_K(eles,E,A,h,s,S)

    nodes = eles+1;
    k_local = (E * A / h) * [1, -1; -1, 1];

    if s==1 then
        // Initialization of stiffness matices of subdomain 1
        K_1 = zeros(nodes, nodes);
        K_s = zeros(nodes-1, nodes-1);
        for i = 1:eles
            K_1(i:i+1, i:i+1) = K_1(i:i+1, i:i+1) + k_local;
        end
        K_s = K_1(2:$,2:$);
    elseif s==S then 
        // Initialization of stiffness matices of subdomain S
        K_s = zeros(nodes, nodes);
        for i = 1 : eles-1
            K_s(i, i) = K_s(i, i) + 2;
            K_s(i, i+1) = K_s(i, i+1) - 1;
            K_s(i+1, i) = K_s(i+1, i) - 1;
        end
        K_s(eles, eles) = K_s(eles, eles) + 1;
        K_s(nodes, nodes) = K_s(nodes, nodes) + 1;
        K_s(1, nodes) = K_s(1, nodes) - 1;
        K_s(nodes, 1) = K_s(nodes, 1) - 1;
        K_s = K_s * (E * A / h);
    elseif (s > 1 && s < S) then 
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
    end

endfunction

//Extract different parts of K local
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

//Create Primal assemly
function AA=build_AA(S)
    AA = zeros((S-1,2*S-2));
    for i = 2:S-1
        AA(i-1,2*i-2) = 1;
        AA(i,2*i-1) = 1;
    end
    AA(1,1)=1;
    AA(S-1,2*S-2)=1;
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
function bbd=build_bbd(S,Fd)
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

/*Global e*/
function ee=build_ee(S,Fd)    
    Rb_conca=build_Rb_conca(S);
    bp_conca=build_bp_conca(S,Fd);
    ee = Rb_conca'*bp_conca;
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


// block diagonal matrix
function A = block_diag(B_list)
    total_rows = 0;
    total_cols = 0;
    for i = 1:length(B_list)
        [rows, cols] = size(B_list(i));
        total_rows = total_rows + rows;
        total_cols = total_cols + cols;
    end

    A = zeros(total_rows, total_cols);

    row_offset = 0;
    col_offset = 0;
    for i = 1:length(B_list)
        [rows, cols] = size(B_list(i));
        A(row_offset+1:row_offset+rows, col_offset+1:col_offset+cols) = B_list(i);
        row_offset = row_offset + rows;
        col_offset = col_offset + cols;
    end
endfunction
