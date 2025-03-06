exec("Tools.sce");

/*********************************
* Section 2 : Domain Decomposition
* Primal Schur Approach - Direct
**********************************/
function Up=Primal_direct(eles,S,E,A,h,Fd)
    Sp_conca=build_Sp_conca(S);
    AA=build_AA(S);
    bp_conca=build_bp_conca(S,Fd);
    SSp = AA*Sp_conca*AA'
    BBp = AA*bp_conca
    Up = SSp\BBp
endfunction

/*********************************
* Section 2 : Domain Decomposition
* Dual Schur Approach - Direct
**********************************/

function Ud=Dual_direct(eles,S,E,A,h,Fd)
    AA_=build_AA_(S);
    Sd_conca=build_Sd_conca(S);
    Rb_conca=build_Rb_conca(S);
    bp_conca=build_bp_conca(S,Fd);
    // Dual formulation of the interface problem 
    SSd = AA_ * Sd_conca * AA_';
    G=build_G(S);
    bbd= AA_*Sd_conca*bp_conca;
    ee=build_ee(S,Fd);
    
    [SSd_r, SSd_c] = size(SSd);
    [G_r, G_c] = size(G);
    [Gt_r, Gt_c] = size(G');
    [bbd_r, bbd_c] = size(bbd);
    [ee_r, ee_c] = size(ee);
    
    dual_M = zeros((SSd_r+Gt_r, SSd_c+G_c));
    dual_rhs = zeros((bbd_r+ee_r,1));
    
    dual_M(1:SSd_r, 1:SSd_c)= SSd;
    dual_M(1:SSd_r, SSd_c+1:$)= G;
    dual_M(SSd_r+1:$, 1:SSd_c)= G';
    dual_rhs(1:bbd_r) = -bbd;
    dual_rhs(bbd_r+1:$) = -ee;
    
    dual_Sol = dual_M\dual_rhs;
    Ud = dual_Sol(bbd_r+1:$);
endfunction
