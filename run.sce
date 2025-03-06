// define total number of elements and number of elements in each subdomain
N_list = [60, 120, 180, 240];
eles_list = [5, 10, 15, 20];

// define length of the bar
L_list = [1.0, 2.0, 3.0, 4.0];

// test values
for N = N_list
    for eles = eles_list
        for L = L_list
            exec("main.sce");
        end
    end
end