function [npcr_avg, uaci_avg] = analyzeKeySensitivity(...
            originalImage, ...
            encrypt_func, keygen_func, dna_decode_func, npcr_uaci_func,...
            initialConditions, rho, sigma, beta, ...
            delta)
% Analyzes key sensitivity by comparing ciphertexts generated with slightly
% different keys using NPCR/UACI.
% Input:
%   originalImage - The image to encrypt (uint8 HxW or HxWx3)
%   encrypt_func - Handle to the encryption function (@encryptImageDNA)
%   keygen_func - Handle to the key generation function (@generateLorenzKeyStream)
%   dna_decode_func - Handle to the DNA decode function (@dna_decode)
%   npcr_uaci_func - Handle to the NPCR/UACI calculation function (@calculateNPCR_UACI)
%   initialConditions, rho, sigma, beta - Original key parameters
%   delta - The small change to apply to one initial condition (e.g., 1e-14)
% Output:
%   npcr_avg - Average NPCR between the two ciphertexts (%)
%   uaci_avg - Average UACI between the two ciphertexts (%)

    fprintf('Running Key Sensitivity Analysis...\n');
    originalSize = size(originalImage);
    rows = originalSize(1);
    cols = originalSize(2);
    numChannels = size(originalImage, 3);
    numPixelsPerChannel = rows * cols;

    % --- 1. Encrypt with Original Key ---
    fprintf('  Encrypting with original key...\n');
    [keyStreamRule_orig, keyStreamDiff_orig] = keygen_func(initialConditions, rho, sigma, beta, numPixelsPerChannel);
    [encryptedDna_orig, ~] = encrypt_func(originalImage, keyStreamRule_orig, keyStreamDiff_orig);
    % Get uint8 version for comparison
    encryptedUint8_orig = getUint8FromDna(encryptedDna_orig, keyStreamRule_orig, dna_decode_func, originalSize);
    fprintf('  Original encryption done.\n');

    % --- 2. Encrypt with Modified Key ---
    initialConditions_modified = initialConditions;
    initialConditions_modified(1) = initialConditions_modified(1) + delta; % Modify first initial condition slightly
    fprintf('  Encrypting with modified key (x0 + %.0e)...\n', delta);
    [keyStreamRule_mod, keyStreamDiff_mod] = keygen_func(initialConditions_modified, rho, sigma, beta, numPixelsPerChannel);
    [encryptedDna_mod, ~] = encrypt_func(originalImage, keyStreamRule_mod, keyStreamDiff_mod);
     % Get uint8 version for comparison
    encryptedUint8_mod = getUint8FromDna(encryptedDna_mod, keyStreamRule_mod, dna_decode_func, originalSize);
    fprintf('  Modified encryption done.\n');

    % --- 3. Compare Ciphertexts using NPCR/UACI ---
    fprintf('  Comparing the two encrypted images...\n');
    [npcr_avg, uaci_avg] = npcr_uaci_func(encryptedUint8_orig, encryptedUint8_mod);

    fprintf('Key sensitivity analysis complete.\n');

end

% --- Helper function for analyzeKeySensitivity ---
function encryptedImageUint8 = getUint8FromDna(encryptedDnaData, keyStreamBinaryRule, dna_decode_func, originalSize)
% Decodes DNA (without inverse diffusion) to get uint8 representation.
    rows = originalSize(1);
    cols = originalSize(2);
    numChannels = originalSize(3);
    numPixelsPerChannel = rows * cols;
    encryptedImageUint8 = zeros(originalSize, 'uint8');

    try
        if numChannels == 1
            decodedBinaryString = dna_decode_func(encryptedDnaData, keyStreamBinaryRule, numPixelsPerChannel);
            decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
            encryptedImageUint8(:) = uint8(bin2dec(decodedBinaryMatrix));
            encryptedImageUint8 = reshape(encryptedImageUint8, rows, cols);
        else % Color
            for k = 1:numChannels
                decodedBinaryString = dna_decode_func(encryptedDnaData{k}, keyStreamBinaryRule, numPixelsPerChannel);
                decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
                channel_uint8 = uint8(bin2dec(decodedBinaryMatrix));
                encryptedImageUint8(:,:,k) = reshape(channel_uint8, rows, cols);
            end
        end
    catch ME
        error('Helper function getUint8FromDna failed: %s', ME.message);
    end
end
