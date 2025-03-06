// define total number of elements and number of elements in each subdomain
N_list = [50, 100];
eles_list = [5, 10, 25];

// define length of the bar
L_list = [1.0, 5.0];

// test values
for N = N_list
    for eles = eles_list
        for L = L_list
            exec("main.sce");
        end
    end
end
