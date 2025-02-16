function cond_vals = femResolution(N,h)
    // Input parameters
    n = N + 1;         // Number of nodes
    E = 1e4;           // Young's Modulus
    A = 1.0;           // Cross-sectional area

    // Local stiffness matrix for one element
    k_local = (E * A / h) * [1, -1; -1, 1];

    // Assemble the global stiffness matrix (K_global)
    K_global = zeros(n, n);
    for i = 1:n-1
        K_global(i:i+1, i:i+1) = K_global(i:i+1, i:i+1) + k_local;
    end

    // Direct Elimination
    K_trim = K_global(2:$, 2:$);

    // Penalty Method
    M = (E * A / h) * 1e4; // Penalty value
    K_penal = K_global;
    K_penal(1,1) = K_penal(1,1) + M;

    // Lagrange Multiplier
    K_extend = zeros(n+1, n+1);
    K_extend(1:n, 1:n) = K_global;
    K_extend(1, n+1) = 1;
    K_extend(n+1, 1) = 1;

    // Compute condition numbers for each method
    cond_trim   = cond(K_trim);
    cond_penal  = cond(K_penal);
    cond_extend = cond(K_extend);

    cond_vals = [cond_trim, cond_penal, cond_extend];
endfunction

// Number of elements
N = 100;
// range of element lengths (h) to test
h_values = linspace(0.001, 10, 100); 

cond_direct   = zeros(1, length(h_values));  // For direct elimination (K_trim)
cond_penalty  = zeros(1, length(h_values));  // For penalty method (K_penal)
cond_lagrange = zeros(1, length(h_values));  // For Lagrange multiplier (K_extend)

// Loop over the different h values and compute condition numbers
for i = 1:length(h_values)
    conds = femResolution(N,h_values(i));
    cond_direct(i)   = conds(1);
    cond_penalty(i)  = conds(2);
    cond_lagrange(i) = conds(3);
end

figure();
plot(h_values, cond_direct, "r");
plot(h_values, cond_penalty, "b");
plot(h_values, cond_lagrange, "g");

xlabel('Element Length (h)');
ylabel('Condition');
title('Variation of Condition with Element Length (h)');
legend('Direct Elimination', 'Penalty Method', 'Lagrange Multiplier');
