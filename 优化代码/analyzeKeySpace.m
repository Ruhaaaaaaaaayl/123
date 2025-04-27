function analyzeKeySpace(initialConditions, parameters)
% Analyzes and prints the theoretical key space based on key components.
% Input:
%   initialConditions - Vector of initial conditions used [x0, y0, z0]
%   parameters - Vector of Lorenz parameters used [rho, sigma, beta] (optional)

    fprintf('Key Components:\n');
    fprintf('  - Lorenz Initial Conditions (x0, y0, z0): [%.15g, %.15g, %.15g]\n', initialConditions(1), initialConditions(2), initialConditions(3));
    if nargin > 1 && ~isempty(parameters)
        fprintf('  - Lorenz Parameters (rho, sigma, beta): [%.4g, %.4g, %.4g]\n', parameters(1), parameters(2), parameters(3));
        fprintf('    (Note: Parameters are often fixed, but can be part of the key)\n');
    end

    fprintf('Theoretical Key Space Discussion:\n');
    fprintf('  - The primary key is typically the set of initial conditions.\n');
    fprintf('  - If using standard IEEE 754 double-precision floating points (64 bits):\n');
    num_initial_conditions = length(initialConditions);
    theoretical_bits_per_condition = 64;
    total_theoretical_bits = num_initial_conditions * theoretical_bits_per_condition;
    fprintf('    - Each initial condition represents %d bits.\n', theoretical_bits_per_condition);
    fprintf('    - Total theoretical bits from initial conditions: %d * %d = %d bits.\n', num_initial_conditions, theoretical_bits_per_condition, total_theoretical_bits);
    fprintf('    - Theoretical key space size: 2^%d (approx 10^%.1f).\n', total_theoretical_bits, total_theoretical_bits * log10(2));
    fprintf('  - Effective Key Space:\n');
    fprintf('    - The *effective* key space might be smaller due to:\n');
    fprintf('      * The precision required for the chaotic system to produce distinct sequences.\n');
    fprintf('      * Potential limitations in the algorithm using the chaotic data.\n');
    fprintf('      * The range of valid initial conditions that maintain chaotic behavior.\n');
    fprintf('    - A common estimate for Lorenz/similar systems suggests an effective precision requirement around 10^-14 to 10^-16, aligning well with double precision.\n');
    fprintf('  - Conclusion: If initial conditions are the key, the space is generally considered very large (> 2^100), sufficient to resist brute-force attacks.\n');

    if nargin > 1 && ~isempty(parameters) && any(parameters ~= [28, 10, 8/3]) % Check if non-standard params used
        fprintf('  - If parameters (rho, sigma, beta) are also part of the key, the key space increases further, but chaotic regions must be considered.\n');
    end
end
