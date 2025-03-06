exec("Section1_FE.sce");
exec("direct_solvers.sce");
exec("iterative_solvers.sce");
exec("Tools.sce");


/*** Input Parameters ***/

// fixed values
E = 1e4;                   // Young'A Modulus
A = 1.0;                   // Surface
Fd = 500;                  // Applied force at the last element (N)


// relative error of two vectors
function err = relative_error(x_ref, x)
    err = norm(x - x_ref) / norm(x_ref);
endfunction

// extract u_ref for interfacial nodes
function u_ref = extract_u_ref(u, S, eles)
    u_ref = zeros(S-1, 1);
    for i = 1:S-1
        u_ref(i) = u(i * eles + 1);
    end
endfunction

// define the tolerance and maximum number of iterations
tol = 1e-6;
max_iter = 1000;

n = N + 1;        // Number of nodes
h = L / N;        // Length of the elements

nodes = eles + 1;   // Number of nodes in each subdomain
S = N / eles;      // Number of subdomains
Inodes = S + 1;        // Number of interfacial nodes
H = L / S;         // Length of the subdomains

K = K_global(E, A, h, n);

/* Executions */ 

/* Section1 */ 
[u, f] = elimination(n, K, Fd);
[u_pen, f_pen] = penalty(E, A, h, K, Fd);
[u_L, f_L] =lagrangian(n, K, Fd);

u_ref = extract_u_ref(u, S, eles);

/* Section2,3 - Direct */ 

start1 = timer();
Up=Primal_direct(eles,S,E,A,h,Fd);
end1 = timer();
t_pd = end1 - start1;

//disp("Up:");
//disp(Up);

start2 = timer();
Ud=Dual_direct(eles,S,E,A,h,Fd);
end2 = timer();
t_dd = end2 - start2;

err_p_direct = relative_error(u_ref, Up);
err_d_direct = relative_error(u_ref, Ud);

//disp("Ud:");
//disp(Ud);

/* Section2 - Primal CG */ 

start3 = timer();
[uub, n_pcg]=Primal_Conjugate_Gradient(eles,S,E,A,h,Fd,max_iter,tol);
end3 = timer();
t_p_cg = end3 - start3;

err_p_cg = relative_error(u_ref, uub);

/* Section3 - Primal CG precond */

start4 = timer();
[u_BDD, n_bdd] = Primal_BDD(eles, S, E, A, h, Fd, max_iter, tol);
end4 = timer();
t_p_bdd = end4 - start4;

err_p_bdd = relative_error(u_ref, u_BDD);

/* Section2 - Dual CG precond */ 

start5 = timer();
[u_b_conca, n_dtefi] = Dual_TEFI(eles,S,E,A,h,Fd,max_iter,tol); 
end5 = timer();
t_d_tefi = end5 - start5;

err_d_tefi = relative_error(u_ref, u_b_conca);

//disp("u_b_conca:");
//disp(u_b_conca);




// save number of iterations till convergence
fd1 = mopen("iter.txt", "at"); 
if (N == 60) & (eles == 5) & (L == 1.0) then
    mfprintf(fd1, "N eles L H h h/H prim_cg prim_bdd dual_tefi\n");
end
mfprintf(fd1, "%d %d %f %f %f %f %d %d %d\n", N, eles, L, H, h, h/H, n_pcg, n_bdd, n_dtefi); 
mclose(fd1);


// save the relative error
fd2 = mopen("error.txt", "at"); 
if (N == 60) & (eles == 5) & (L == 1.0) then
    mfprintf(fd2, "N eles L H h h/H prim_dir dual_dir prim_cg prim_bdd dual_tefi\n");
end
mfprintf(fd2, "%d %d %f %f %f %f %e %e %e %e %e\n", N, eles, L, H, h, h/H, err_p_direct, err_d_direct, err_p_cg, err_p_bdd, err_d_tefi); 
mclose(fd2);


// save the time
fd3 = mopen("time.txt", "at");
if (N == 60) & (eles == 5) & (L == 1.0) then
    mfprintf(fd3, "N eles L H h h/H prim_dir dual_dir prim_cg prim_bdd dual_tefi\n");
end
mfprintf(fd3, "%d %d %f %f %f %f %f %f %f %f %f\n", N, eles, L, H, h, h/H, t_pd, t_dd, t_p_cg, t_p_bdd, t_d_tefi);
mclose(fd3);


// plot results

x_nodes = (1:n) * h;
x_interface = (1:S-1) * H;

y = Fd / (E * A) * x_nodes; // analytical solution

figure();
plot(x_nodes, y, '-k');  
plot(x_nodes, u, '-b');   
plot(x_nodes, u_pen, '-r'); 
plot(x_nodes, u_L, '-g');  
plot(x_interface, Up, '-c');  
plot(x_interface, Ud, '-m');  
plot(x_interface, uub, '-y'); 
plot(x_interface, u_BDD, '-', color=[1, 0.647, 0]); 
plot(x_interface, u_b_conca, '-', color=[0.5, 0, 0.5]); 

legend(["Analytical solution", "Elimination", "Penalty", "Lagrangian", "Primal direct", "Dual direct", "Primal CG", "Primal BDD", "Dual TEFI"], "in_upper_left");

xlabel("x");
ylabel("u");
title("N=" + string(N) + ", eles=" + string(eles) + ", L=" + string(L));
