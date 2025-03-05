
//Create stiffness matrix
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
function Sp=build_Sp(K, s, S)
    [Kii, Kib, Kbi, Kbb] = extract_K(K, s, S);
    Sp = Kbb - Kbi * inv(Kii) * Kib;
endfunction

/* Concatenated Sp */
function Sp_conca=build_Sp_conca(K, S)
    Sp_conca = zeros((2*S-2,2*S-2));
    Sp_1=build_Sp(K,1,S);
    Sp_S=build_Sp(K,S,S);
    
    if S>2 then
        Sp_s=build_Sp(K,2,S);
        for i = 1:S-2
            Sp_conca(2*i:2*i+1,2*i:2*i+1)= Sp_s;
        end
    end
    Sp_conca(1,1)=Sp_1;
    Sp_conca(2*S-2,2*S-2)=Sp_S;
endfunction

/* Local Sd */
function Sd=build_Sd(K, s, S)
    Sp=build_Sp(K, s, S);
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

//Conjugate gradient
function [uub, nb_iter, rel_err]=Conjugate_Gradient(eles,S,E,A,h,Fd,m,tol)
    //preliminary 
    [K_1_strim K_s K_S]=build_K(eles,E,A,h);
    [k1_r,k1_c] = size(K_1_strim);
    [ks_r,ks_c] = size(K_s);
    [kS_r,kS_c] = size(K_S);
    AA = build_AA(S);
    
    /*** Initialization Step ***/
    //ub0
    uub = zeros(S-1,1);
    //Compure concatenated residual rb = zeros(2*S-2,1);
    rb_concatanated = [];
    //Compute local concatenated vector
    ub_concatanated = AA' * uub;
    
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
        ub0_s = ub_concatanated(2*(s-1):2*(s-1)+1);
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
    ddb=rrb;
    
    //Check onvergence 
    norm_r0 = norm(rrb);  
    if norm_r0 == 0 then
       nb_iter = 0;
       rel_err = 0;
       disp("Initial residual is zero. Already converged.");
       return
    end
    
    /*** Iterations ***/ 
    for k=1:m
        // Compute local vector
        db_concatanated = AA' * ddb ;
        local_matvec_concatanated = [];
        
        //s=1
        [Kii, Kib, Kbi, Kbb] = extract_K(K_1_strim, 1, S);
        db_1 = db_concatanated(1);
        di_1 = -inv(Kii)*Kib*db_1;
        res_1 = K_1_strim*[di_1; db_1] ;
        local_matvec_concatanated = [local_matvec_concatanated; res_1($)];
        
        //s=2...S-1
        for s=2:S-1
            [Kii, Kib, Kbi, Kbb] = extract_K(K_s, s, S);
            db_s = db_concatanated(2*(s-1):2*(s-1)+1);
            di_s = -inv(Kii)*Kib*db_s;
            res_s = K_s*[di_s; db_s];
            local_matvec_concatanated = [local_matvec_concatanated; res_s($-1:$)];
        end
     
        //s=S
        [Kii, Kib, Kbi, Kbb] = extract_K(K_S, S, S);
        db_S = db_concatanated($);
        di_S = -inv(Kii)*Kib*db_S;
        res_S = K_S*[di_S; db_S] ;
        local_matvec_concatanated = [local_matvec_concatanated; res_S($)];
        
        //Compute global matrix vector product 
        SSdd = AA*local_matvec_concatanated;
        
        //Compute optimal step
        alpha = (rrb'*ddb)/(ddb'*SSdd);
        
        //Compute iterate
        uub = uub+alpha*ddb;
        
        //Compute residual
        rrb = rrb - alpha*SSdd;
                            
        //Compute orthogonalization parameter
        beta1 = - (rrb'*SSdd)/(ddb'*SSdd);
        
        //Update search direction
        ddb = rrb + beta1*ddb ;
        
        //Check onvergence 
        rel_error = norm(rrb) / norm_r0; 
        if rel_error < tol then
           nb_iter = k;
           rel_err = rel_error;
           disp("Converged at iteration " + string(k));
           break;
        end
    end
    
    nb_iter = m;
    rel_err = norm(rrb) / norm_r0;
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

// preconditioned primal CG
// primal schur BDD
function [u, nb_iter, rel_err] = primal_schur_BDD(eles, S, E, A, h, Fd, m, tol)
    
    [K_1_strim, K_s, K_S] = build_K(eles, E, A, h);
    AA = build_AA(S);
    disp("AA: ");
    disp(AA);

    G = build_G(S); 
    disp("G: ");
    disp(G);

    /*** initialization ***/
    bp_conca = [];
    Sp_list = list();
    for s = 1:S
        if s == 1 then
            [Kii, Kib, Kbi, Kbb] = extract_K(K_1_strim, s, S);
            fi = zeros(size(Kii, 1), 1);
            fb = 0;
        elseif s == S then
            [Kii, Kib, Kbi, Kbb] = extract_K(K_S, s, S);
            fi = zeros(size(Kii, 1), 1);
            fb = Fd;
        else
            [Kii, Kib, Kbi, Kbb] = extract_K(K_s, s, S);
            fi = zeros(size(Kii, 1), 1);
            fb = zeros(2, 1);
        end
        Sp = Kbb - Kbi * inv(Kii) * Kib;
        Sp_list(s) = Sp;
        bp_conca = [bp_conca; fb];
    end

    Sp_conca = block_diag(Sp_list);
    disp("Sp_conca: ");
    disp(Sp_conca);
    
    Sp_a = AA' * Sp_conca * AA;
    Sd = pinv(Sp_a);
    P = eye(S-1, S-1) - G * inv(G' * Sp_a * G) * G' * Sp_a;

    u = G * inv(G' * Sp_a * G) * G' * bp_conca; 
    r = P' * bp_conca;
    z = Sd * r;
    d = [];
    d(1) = z
    p = [];
    
    // 
    norm_r0 = norm(r);
    if norm_r0 == 0 then
        nb_iter = 0;
        rel_err = 0;
        disp("Initial residual is zero. Already converged.");
        return;
    end

    /*** iteration ***/
    beta = 0;
    for k = 1:m+1
        p(k) = P' * Sp_a * d(k);
        alpha = (r' * d(k)) / (d(k)' * p(k));
        u = u + alpha * d(k);
        r = r - alpha * p(k);
        z = Sd_conca * r;

        for j = 1:k
            beta =  beta - (z' *  p(j)) / (d(j)' * p(j));
        end
        d(k + 1) = z + beta * d(k);

        rel_error = norm(r) / norm_r0;
        if rel_error < tol then
            nb_iter = k;
            rel_err = rel_error;
            disp("Converged at iteration " + string(k));
            return;
        end
    end

    //
    nb_iter = m;
    rel_err = norm(r) / norm_r0;
    disp("Max iteration reached. Not converged.");
endfunction