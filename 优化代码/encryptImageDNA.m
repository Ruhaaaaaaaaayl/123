function [encryptedDnaData, originalSize] = encryptImageDNA(originalImage, keyStreamBinaryRule, keyStreamByteDiffusion)
% Encrypts Grayscale or Color Image: Diffusion -> DNA Encoding
% Input:
%   originalImage - Original image matrix (uint8 HxW or HxWx3)
%   keyStreamBinaryRule - Key for DNA rules (char 1 x (N_per_channel*4))
%   keyStreamByteDiffusion - Key for diffusion (uint8 1 x N_per_channel)
% Output:
%   encryptedDnaData - Encrypted DNA sequence(s).
%                      Char array (1x(N*4)) for grayscale.
%                      Cell array {1x3} of char arrays for color.
%   originalSize - Original image size [rows, cols, channels]

    if ~isa(originalImage, 'uint8')
        error('Input image must be uint8 type.');
    end
    originalSize = size(originalImage);
    rows = originalSize(1);
    cols = originalSize(2);
    numChannels = size(originalImage, 3);
    numPixelsPerChannel = rows * cols;

    % --- Check Key Stream Lengths against per-channel size ---
    if length(keyStreamBinaryRule) ~= numPixelsPerChannel * 4
        error('Rule key stream length (%d) incorrect for %d pixels (expected %d).', length(keyStreamBinaryRule), numPixelsPerChannel, numPixelsPerChannel * 4);
    end
    if length(keyStreamByteDiffusion) ~= numPixelsPerChannel
         error('Diffusion key stream length (%d) incorrect for %d pixels.', length(keyStreamByteDiffusion), numPixelsPerChannel);
    end

    % --- Process each channel ---
    if numChannels == 1
        % --- Grayscale Processing ---
        img_vector = originalImage(:)'; % Flatten to 1 x N

        % 1. Diffusion Layer (CFB-like XOR)
        diffused_vector = zeros(1, numPixelsPerChannel, 'uint8');
        diffused_vector(1) = bitxor(img_vector(1), keyStreamByteDiffusion(1)); % C1 = P1 XOR K1
        for i = 2:numPixelsPerChannel
            feedback = diffused_vector(i-1); % Use previous *ciphertext* byte
            key_feedback_xor = bitxor(keyStreamByteDiffusion(i), feedback);
            diffused_vector(i) = bitxor(img_vector(i), key_feedback_xor); % Ci = Pi XOR (Ki XOR Ci-1)
        end

        % 2. Substitution Layer (DNA Encoding)
        encryptedDnaData = dna_encode(diffused_vector, keyStreamBinaryRule, numPixelsPerChannel);

    elseif numChannels == 3
        % --- Color Processing (Channel by Channel) ---
        encryptedDnaData = cell(1, 3); % Initialize cell array for results
        for k = 1:3
            img_vector = originalImage(:,:,k); % Extract channel
            img_vector = img_vector(:)'; % Flatten to 1 x N

            % 1. Diffusion Layer (using the same key stream for each channel)
            diffused_vector = zeros(1, numPixelsPerChannel, 'uint8');
            diffused_vector(1) = bitxor(img_vector(1), keyStreamByteDiffusion(1));
            for i = 2:numPixelsPerChannel
                feedback = diffused_vector(i-1);
                key_feedback_xor = bitxor(keyStreamByteDiffusion(i), feedback);
                diffused_vector(i) = bitxor(img_vector(i), key_feedback_xor);
            end

            % 2. Substitution Layer (DNA Encoding, same rule key for each channel)
            encryptedDnaData{k} = dna_encode(diffused_vector, keyStreamBinaryRule, numPixelsPerChannel);
            % fprintf('Channel %d encrypted.\n', k); % Optional progress
        end
    else
        error('Unsupported number of image channels: %d', numChannels);
    end
end
