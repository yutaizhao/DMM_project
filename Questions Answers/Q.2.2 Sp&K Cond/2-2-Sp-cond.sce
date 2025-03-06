function sp_cond = spResolution(N, H, h)

    E = 1e4;              // Young's Modulus
    A = 1.0;              // Cross-sectional area
    eles = round(H/h);            // Number of elements per subdomain
    nodes = eles + 1;     // Number of nodes per subdomain
    S = N / eles;         // Number of subdomains (assumed integer)

    k_local = (E * A / h) * [1, -1; -1, 1];

    K_1 = zeros(nodes, nodes);
    for i = 1:eles
        K_1(i:i+1, i:i+1) = K_1(i:i+1, i:i+1) + k_local;
    end
    K_1_strim = K_1(2:$, 2:$);

    K_s = zeros(nodes, nodes);
    for i = 1:eles-1
        K_s(i, i) = K_s(i, i) + 2;
    end
    for i = 1:eles-2
        K_s(i, i+1) = K_s(i, i+1) - 1;
        K_s(i+1, i) = K_s(i+1, i) - 1;
    end
    K_s(1, eles)       = K_s(1, eles) - 1;
    K_s(eles, 1)       = K_s(eles, 1) - 1;
    K_s(nodes, eles-1) = K_s(nodes, eles-1) - 1;
    K_s(eles-1, nodes) = K_s(eles-1, nodes) - 1;
    K_s(eles, eles)    = K_s(eles, eles) + 1;
    K_s(nodes, nodes)  = K_s(nodes, nodes) + 1;
    K_s = K_s * (E * A / h);

    K_S = zeros(nodes, nodes);
    for i = 1:eles-1
        K_S(i, i)     = K_S(i, i) + 2;
        K_S(i, i+1)   = K_S(i, i+1) - 1;
        K_S(i+1, i)   = K_S(i+1, i) - 1;
    end
    K_S(eles, eles)   = K_S(eles, eles) + 1;
    K_S(nodes, nodes) = K_S(nodes, nodes) + 1;
    K_S(1, nodes)     = K_S(1, nodes) - 1;
    K_S(nodes, 1)     = K_S(nodes, 1) - 1;
    K_S = K_S * (E * A / h);

    [k1s_r, k1s_c] = size(K_1_strim);
    [ks_r, ks_c]   = size(K_s);
    [kS_r, kS_c]   = size(K_S);

    // Primal Schur complement for subdomain 1
    Sp_1 = K_1_strim(k1s_r, k1s_c) - K_1_strim(k1s_r, 1:k1s_r-1) * inv(K_1_strim(1:k1s_r-1, 1:k1s_r-1)) * K_1_strim(1:k1s_r-1, k1s_c);

    // Primal Schur complement for an interior subdomain (subdomain s)
    Sp_s = K_s(ks_r-1:ks_r, ks_c-1:ks_c) - K_s(ks_r-1:ks_r, 1:ks_r-2) * inv(K_s(1:ks_r-2, 1:ks_r-2)) * K_s(1:ks_r-2, ks_c-1:ks_c);

    // Primal Schur complement for the last subdomain (subdomain S)
    Sp_S = K_S(kS_r, kS_c) - K_S(kS_r, 1:kS_r-1) * inv(K_S(1:kS_r-1, 1:kS_r-1)) * K_S(1:kS_r-1, kS_c);

    Sp = zeros((2*S - 2, 2*S - 2));
    for i = 1:S-2
        Sp(2*i:2*i+1, 2*i:2*i+1) = Sp_s;
    end
    // Place the first and last subdomain contributions
    Sp(1,1) = Sp_1;
    Sp(2*S-2, 2*S-2) = Sp_S;

    sp_cond = cond(Sp);
endfunction

N=100;

// Define a range of element lengths (h) and range of domain lengths (H)
h_values = linspace(0.001, 10, 100);
H_values = linspace(0.01, 100, 100);

x_values = h_values ./ H_values;

// Initialize a vector to store the condition numbers of Sp
cond_Sp = zeros(length(h_values));

// Loop over the different h values and compute cond(Sp)
for i = 1:length(h_values)
    cond_Sp(i) = spResolution(N, H_values(i), h_values(i));
end

// Plot the variation of the condition number of Sp versus h/H
figure();
scatter(x_values, cond_Sp);  // Magenta line with circle markers
xlabel("h/H");
ylabel("Condition of Sp");
title("Variation of Sp Condition vs. h/H");
legend("Sp Condition", "location", "best");
