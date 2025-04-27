function decryptedImage = decryptImageDNA(encryptedDnaData, originalSize, keyStreamBinaryRule, keyStreamByteDiffusion)
% Decrypts Grayscale or Color Image: DNA Decoding -> Inverse Diffusion
% Input:
%   encryptedDnaData - Encrypted DNA sequence(s).
%                      Char array (1x(N*4)) for grayscale.
%                      Cell array {1x3} of char arrays for color.
%   originalSize - Original image size [rows, cols, channels]
%   keyStreamBinaryRule - Key for DNA rules (char 1 x (N_per_channel*4))
%   keyStreamByteDiffusion - Key for diffusion (uint8 1 x N_per_channel)
% Output:
%   decryptedImage - Decrypted image matrix (uint8 HxW or HxWx3)

    rows = originalSize(1);
    cols = originalSize(2);
    numChannels = originalSize(3); % Get channels from originalSize
    numPixelsPerChannel = rows * cols;

    % --- Check Key Stream Lengths ---
     if length(keyStreamBinaryRule) ~= numPixelsPerChannel * 4
        error('Rule key stream length (%d) incorrect for %d pixels (expected %d).', length(keyStreamBinaryRule), numPixelsPerChannel, numPixelsPerChannel * 4);
    end
    if length(keyStreamByteDiffusion) ~= numPixelsPerChannel
         error('Diffusion key stream length (%d) incorrect for %d pixels.', length(keyStreamByteDiffusion), numPixelsPerChannel);
    end

    % --- Initialize output ---
    decryptedImage = zeros(originalSize, 'uint8');

    % --- Process based on number of channels ---
    if numChannels == 1
        % --- Grayscale Decryption ---
        if ~ischar(encryptedDnaData)
            error('Expected char array for grayscale DNA data, received %s.', class(encryptedDnaData));
        end
         if length(encryptedDnaData) ~= numPixelsPerChannel * 4
             error('Encrypted grayscale DNA sequence length (%d) is incorrect for %d pixels (expected %d).', length(encryptedDnaData), numPixelsPerChannel, numPixelsPerChannel * 4);
         end

        % 1. Inverse Substitution (DNA Decoding)
        decodedBinaryString = dna_decode(encryptedDnaData, keyStreamBinaryRule, numPixelsPerChannel);
        try
            decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
            diffused_vector = uint8(bin2dec(decodedBinaryMatrix))'; % This is the state AFTER diffusion
        catch ME
            error('Error converting decoded binary to uint8 (grayscale): %s', ME.message);
        end

        % 2. Inverse Diffusion Layer
        % P1 = C1 XOR K1
        % Pi = Ci XOR (Ki XOR Ci-1) (where C is the diffused_vector)
        decrypted_vector = zeros(1, numPixelsPerChannel, 'uint8');
        decrypted_vector(1) = bitxor(diffused_vector(1), keyStreamByteDiffusion(1));
        for i = 2:numPixelsPerChannel
            feedback = diffused_vector(i-1); % Use previous *diffused* byte (Ci-1)
            key_feedback_xor = bitxor(keyStreamByteDiffusion(i), feedback);
            decrypted_vector(i) = bitxor(diffused_vector(i), key_feedback_xor);
        end

        % 3. Reshape
        try
            decryptedImage = reshape(decrypted_vector, rows, cols);
        catch ME
             error('Error reshaping decrypted grayscale vector: %s. Size mismatch?', ME.message);
        end

    elseif numChannels == 3
        % --- Color Decryption (Channel by Channel) ---
        if ~iscell(encryptedDnaData) || numel(encryptedDnaData) ~= 3
             error('Expected cell array of size 3 for color DNA data, received %s.', class(encryptedDnaData));
        end

        for k = 1:3
            dna_seq = encryptedDnaData{k};
            if ~ischar(dna_seq) || length(dna_seq) ~= numPixelsPerChannel * 4
                error('Encrypted color DNA sequence for channel %d has incorrect type or length.', k);
            end

            % 1. Inverse Substitution (DNA Decoding)
            decodedBinaryString = dna_decode(dna_seq, keyStreamBinaryRule, numPixelsPerChannel);
             try
                decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
                diffused_vector = uint8(bin2dec(decodedBinaryMatrix))'; % State AFTER diffusion for this channel
            catch ME
                error('Error converting decoded binary to uint8 (channel %d): %s', k, ME.message);
            end

            % 2. Inverse Diffusion Layer (using same keys)
            decrypted_vector = zeros(1, numPixelsPerChannel, 'uint8');
            decrypted_vector(1) = bitxor(diffused_vector(1), keyStreamByteDiffusion(1));
             for i = 2:numPixelsPerChannel
                feedback = diffused_vector(i-1);
                key_feedback_xor = bitxor(keyStreamByteDiffusion(i), feedback);
                decrypted_vector(i) = bitxor(diffused_vector(i), key_feedback_xor);
             end

             % 3. Reshape and place in output image
            try
                decryptedImage(:,:,k) = reshape(decrypted_vector, rows, cols);
            catch ME
                 error('Error reshaping decrypted vector for channel %d: %s', k, ME.message);
            end
            % fprintf('Channel %d decrypted.\n', k); % Optional progress
        end
    else
        error('Unsupported number of image channels in originalSize: %d', numChannels);
    end
end

