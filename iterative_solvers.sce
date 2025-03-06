exec("Tools.sce");

/*Please note that we are in the case of homogeneuos structure so m and Q=I*/

    /*********************************
    * Section 2 : Domain Decomposition
    * Primal Schur Approach - CG
    **********************************/

    function [uub, n_iter] =Primal_Conjugate_Gradient(eles,S,E,A,h,Fd,m,tol)

        //preliminary 

        AA = build_AA(S);

        /*** Initialization Step ***/

        //ub0
        uub = zeros(S-1,1);
        //Compure concatenated residual rb = zeros(2*S-2,1);
        rb_concatanated = [];
        //Compute local concatenated vector
        ub_concatanated = AA' * uub;

        //s=1
        K_1_strim=build_K(eles,E,A,h,1,S);
        [k1_r,k1_c] = size(K_1_strim);
        [Kii, Kib, Kbi, Kbb] = extract_K(K_1_strim, 1, S);
        fi_1=zeros(k1_r-1,1);
        fb_1=0;
        ub0_1 = ub_concatanated(1,1);
        ui0_1 = inv(Kii)*(fi_1-Kib*ub0_1);
        res_1 = K_1_strim*[ui0_1; ub0_1]-[fi_1; fb_1] ;
        rb_concatanated = [rb_concatanated; -res_1($)];

        //s=2...S-1
        K_s=build_K(eles,E,A,h,2,S);
        [ks_r,ks_c] = size(K_s);
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
        K_S=build_K(eles,E,A,h,S,S);
        [kS_r,kS_c] = size(K_S);
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
            n_iter = 0;
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
                n_iter = k;
                disp("Converged at iteration " + string(k));
                break;
            end
        end

    endfunction


    /*********************************
    * Section 2 : Domain Decomposition
    * Primal BDD - Preconditioned CG
    **********************************/

    function [u, n_iter]= Primal_BDD(eles, S, E, A, h, Fd, m, tol)

        //Preliminary 

        /*P*/
        AA=build_AA(S);
        Sp_conca=build_Sp_conca(S);
        SSp = AA*Sp_conca*AA'
        G = build_G(S);
        P = eye(S-1,S-1) - G*inv(G'*SSp*G)*G'*SSp;

        /*u0*/
        bp_conca=build_bp_conca(S,Fd);
        BBp = AA*bp_conca;
        u0 = G*inv(G'*SSp*G)*G'*BBp;  

        /*** Initialization Step ***/

        //u(0)
        u = u0;

        //Compute local concatenated u
        u_concatanated = AA' * u;

        /** Compute r0 globally **/
        rr = P'*BBp;

        //Check onvergence 
        norm_r0 = norm(rr);  
        if norm_r0 < tol then
            n_iter = 0;
            disp("BDD : Converged at iteration 0");
            /* Post processing */
            return;
        end

        /** Compute z0,d0 **/
        SSp_inv=build_SSp_inv(S);
        zz = SSp_inv*rr;
        dd = zz;

        p_list=[];
        dd_list=[dd];

        /* For loop */

        for i=1:m 

            // Compute local p
            d_concatanated = AA' * dd;
            local_Spxd_concatanated = [];

            Sp = build_Sp(1,S);
            d = d_concatanated(1);
            local_Spxd_concatanated = [local_Spxd_concatanated, Sp*d];
            for s=2:S-1
                Sp = build_Sp(s,S);
                d = d_concatanated(2*(s-1):2*(s-1)+1);
                local_Spxd_concatanated = [local_Spxd_concatanated; Sp*d];
            end 
            Sp = build_Sd(S,S);
            d = d_concatanated($);
            local_Spxd_concatanated = [local_Spxd_concatanated; Sp*d];
            //END of computation p

            p = P'*AA*local_Spxd_concatanated;
            p_list=[p_list;p];

            // Update CG iter
            a = (rr'*dd)/(dd'*p);
            u = u + a*dd;

            rr = rr -a*p;
            zz = SSp_inv*rr;

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
            rel_error = norm(rr) / norm_r0; 
            if rel_error < tol then
                n_iter = i;
                disp("Converged at iteration " + string(i));
                break;
            end

        end

    endfunction



    /*********************************
    * Section 2 : Domain Decomposition
    * Dual TEFI - Preconditioned CG
    **********************************/
    function [u_b_conca, n_iter] = Dual_TEFI(eles, S, E, A, h, Fd, m, tol)  

        //Preliminary 

        /*P*/
        G = build_G(S);
        P = eye(S-1,S-1) - G*inv(G'*G)*G';

        /*lambda0*/
        ee=build_ee(S,Fd);
        lambda0 =  - G*inv(G'*G)*ee; 

        AA_=build_AA_(S);
        SSd_inv=build_SSd_inv(S);

        /*** Initialization Step ***/

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
        if norm_r0 < tol then
            n_iter = 0;
            disp("TEFI : Converged at iteration 0");
            /* Post processing */
            Sd_conca=build_Sd_conca(S);
            SSd = AA_*Sd_conca*AA_';
            bbd=build_bbd(S);
            bp_conca=build_bp_conca(S,Fd);
            Rb_conca=build_Rb_conca(S);
            alpha_b = inv(G'*G)*G'*(-bbd-SSd*Lambda);
            u_b_conca = Sd_conca*(bp_conca+AA_'*Lambda)+Rb_conca*alpha_b;
            return;
        end

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
            a = (rr'*dd)/(dd'*p);
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
            rel_error = norm(rr) /norm_r0; 
            if rel_error < tol then
                n_iter = i;
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
        u_b_conca = Sd_conca*(bp_conca+AA_'*Lambda)+Rb_conca*alpha_b;

    endfunction


