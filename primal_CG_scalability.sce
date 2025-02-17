/*********************************
* Section 2 : Domain Decomposition
* Primal Schur Approach - CG
**********************************/


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



function scalability()
    
    E = 1e4;
    A = 1.0;
    L = 1.0; 
    Fd = 500;
    maxIter = 100;
    tol = 0.0001;

    // h/H = const
    // N = S * eles
    // => h/H = S/N = 1/eles

    eles = 5;
    S_list = [2, 4, 8, 16];
    nbIter_list = zeros(size(S_list));
    
    for i = 1:length(S_list)
        S = S_list(i);
        N = S * eles;
        h = L / N;
        
        [Up, nb_iter, rel_err] = Conjugate_Gradient(eles, S, E, A, h, Fd, maxIter, tol);
        
        nbIter_list(i) = nb_iter;
    end
    
    figure;
    plot(S_list, nbIter_list, 'o-',);
    xlabel('Number of subdomains S');
    ylabel('Number of CG iterations');
    title('Primal Schur CG: Scalability');

end

scalability();



