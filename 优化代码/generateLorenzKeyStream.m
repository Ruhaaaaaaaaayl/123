function [keyStreamBinaryRule, keyStreamByteDiffusion] = generateLorenzKeyStream(initialConditions, rho, sigma, beta, numPixels)
% Uses Lorenz system to generate chaotic key streams for DNA rule selection
% and byte-level diffusion.
% Input:
%   initialConditions - Lorenz initial conditions [x0, y0, z0] (Key!)
%   rho, sigma, beta - Lorenz parameters (can be part of key)
%   numPixels - Number of pixels PER CHANNEL (H * W)
% Output:
%   keyStreamBinaryRule - Binary key stream for DNA rules (char 1 x (numPixels * 4))
%   keyStreamByteDiffusion - Key stream for byte diffusion (uint8 1 x numPixels)

    % --- Parameters ---
    h = 0.001;         % Step size
    transient_steps = 1000; % Initial steps to discard

    % Points needed: numPixels for rules (x), numPixels for diffusion (z)
    required_points = numPixels;
    num_steps = transient_steps + required_points + 100; % Add buffer

    % --- Solve Lorenz System (Euler method) ---
    x = zeros(1, num_steps);
    y = zeros(1, num_steps);
    z = zeros(1, num_steps);
    x(1) = initialConditions(1);
    y(1) = initialConditions(2);
    z(1) = initialConditions(3);

    if any(abs(initialConditions) < 1e-6)
       warning('Initial conditions are very close to zero, Lorenz system might stabilize or behave predictably.');
    end
    if any(isnan(initialConditions)) || any(isinf(initialConditions))
        error('Initial conditions cannot be NaN or Inf.');
    end


    for i = 1:num_steps-1
        dx = sigma * (y(i) - x(i));
        dy = x(i) * (rho - z(i)) - y(i);
        dz = x(i) * y(i) - beta * z(i);

        x(i+1) = x(i) + h * dx;
        y(i+1) = y(i) + h * dy;
        z(i+1) = z(i) + h * dz;

        if ~all(isfinite([x(i+1), y(i+1), z(i+1)])) || any(abs([x(i+1), y(i+1), z(i+1)]) > 1e6)
             error('Lorenz system diverged or produced non-finite values. Check parameters (rho=%.2f, sigma=%.2f, beta=%.2f) and initial conditions [%.2e, %.2e, %.2e].', rho, sigma, beta, initialConditions(1), initialConditions(2), initialConditions(3));
        end
    end

    % --- Discard transient part ---
    x_steady = x(transient_steps+1 : end);
    z_steady = z(transient_steps+1 : end);

    available_points = length(x_steady);
    if available_points < required_points
        error('Generated Lorenz steady sequence length (%d) is insufficient, need %d points. Increase num_steps.', available_points, required_points);
    end

    % --- Generate Key for DNA Rules (using x) ---
    % Scale, take mod 16, convert to 4-bit binary string
    key_rule_int = mod(floor(abs(x_steady(1:numPixels)) * 1e8), 16);
    keyStreamBinaryRule = dec2bin(key_rule_int, 4); % N x 4 char array
    keyStreamBinaryRule = keyStreamBinaryRule';     % 4 x N
    keyStreamBinaryRule = keyStreamBinaryRule(:)';  % 1 x (N*4) char row vector

    % --- Generate Key for Byte Diffusion (using z) ---
    % Scale, take mod 256, convert to uint8
    keyStreamByteDiffusion = uint8(mod(floor(abs(z_steady(1:numPixels)) * 1e8), 256));
    keyStreamByteDiffusion = keyStreamByteDiffusion(:)'; % Ensure 1 x N uint8 row vector

    % --- Final Length Check ---
    if length(keyStreamBinaryRule) ~= numPixels * 4
        error('Internal error: Rule key stream length mismatch (%d vs %d).', length(keyStreamBinaryRule), numPixels*4);
    end
    if length(keyStreamByteDiffusion) ~= numPixels
        error('Internal error: Diffusion key stream length mismatch (%d vs %d).', length(keyStreamByteDiffusion), numPixels);
    end
end
