/*********************************
* Section 3 : Domain Decomposition
* Mixed approach - LATIN Method
**********************************/

/*** Input Parameters ***/
L = 1.0;          // Total length of the domain
E = 1e4;          // Young's Modulus
A = 1.0;          // Cross-sectional area
Fd = 500;         // Applied force at the right end (N)
N = 20;           // Total number of elements in the domain
h = L / N;        // Length of each element

eles_sub = 5;         // Number of elements per subdomain (can be adjusted)
S = N / eles_sub;     // Number of subdomains (assume N is divisible by eles_sub)
Inodes = S + 1;       // Number of interface nodes

k0 = E * A / h;       // Element stiffness coefficient


function K = build_K(eles, k0)
    nodes = eles + 1;
    k_elem = k0 * [1, -1; -1, 1];
    K = zeros(nodes, nodes);
    for i = 1:eles
        K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + k_elem;
    end
endfunction

function k = build_k(eles, inter_k0, s, S)
    nodes = eles + 1;
    k = zeros(nodes, nodes);
    k(1,1) = inter_k0(s);
    k(nodes, nodes) = inter_k0(s+1);
endfunction


function local_list = LATIN(eles, S, E, A, h, Fd, m, tol)
    
    // Define the stiffness at the interfaces
    inter_k0 = zeros(1, S+1);
    inter_k0(1) = 1000 * k0; 
    for s = 2:S
        inter_k0(s) = k0 / (s-1);
    end
    inter_k0(S+1) = 0; 
    
    inv_sum_list = list();
    for s = 1:S
        K_sub = build_K(eles, k0);
        k_sub = build_k(eles, inter_k0, s, S);
        A_sub = K_sub + k_sub;  
        inv_A_sub = inv(A_sub); 
        inv_sum_list($+1) = inv_A_sub;
    end
    
    // local_list rows are [W_hat, F_hat]
    // There are S+1 interface nodes.
    local_list = zeros(S+1, 2);
    // Boundary conditions: left end displacement = 0, right end force = Fd.
    local_list(1, :) = [0, 1000 * Fd/k0];
    local_list(S+1, :) = [0, Fd];
    
    for iter = 1:m
       
        linear_list = zeros(2*S, 2);
        
        // Linear stage
        for s = 1:S
            // Only the first and last nodes are driven by interface data.
            b = zeros(eles+1, 1);
            b(1) = - local_list(s, 2) + inter_k0(s) * local_list(s, 1);
            b($) = local_list(s+1, 2) + inter_k0(s+1) * local_list(s+1, 1);
            
            U = inv_sum_list(s) * b;
            // Extract the interface displacements (first and last node)
            W = [ U(1); U($) ];
            F = [ - local_list(s, 2) + inter_k0(s)*(local_list(s, 1) - W(1)); local_list(s+1, 2) + inter_k0(s+1)*(local_list(s+1, 1) - W(2)) ];
            linear_list(2*s-1, :) = [W(1), F(1)];
            linear_list(2*s, :)   = [W(2), F(2)];
        end
        
        // Local stage
        new_local_list = zeros(size(local_list));
        new_local_list(1, :) = [0, linear_list(1, 2)];  // Left boundary: displacement fixed at 0
        for s = 2:S
            left = 2*s - 2;
            right = 2*s - 1;
            W_hat = 0.5*(linear_list(left, 1) + linear_list(right, 1)) - 0.5*(1 / inter_k0(s))*(linear_list(left, 2) + linear_list(right, 2));
            F_hat = linear_list(left, 2) + inter_k0(s)*(W_hat - linear_list(left, 1));
            new_local_list(s, :) = [W_hat, F_hat];
        end
        new_local_list(S+1, :) = [linear_list($, 1), Fd]; // Right boundary: force is Fd
        
        // Check convergence (using norm instead of computing eta - we don't know how to do) 
        if norm(new_local_list - local_list) < tol then
            local_list = new_local_list;
            disp("Converged at iteration " + string(iter));
            break;
        end
        
        local_list = new_local_list;
    end
endfunction

// Call the LATIN method without relaxation
local_list = LATIN(eles_sub, S, E, A, h, Fd, 5000, 0.0001);
u_b = local_list(:, 1);
disp("Interface displacements by LATIN :");
disp(u_b);
