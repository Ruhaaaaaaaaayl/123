function [decryptedNoisyImage, psnr_val, mse_val] = analyzeNoiseAttack(...
            originalImage, encryptedImageUint8, ...
            keyStreamBinaryRule, keyStreamByteDiffusion, ...
            noiseType, noiseParam, ...
            dna_encode_func, decryptImageDNA_func)
% Simulates a noise attack on the encrypted image and decrypts it.
% Input:
%   originalImage - The original clean image (for PSNR comparison)
%   encryptedImageUint8 - The encrypted image (uint8 HxW or HxWx3) to attack
%   keyStreamBinaryRule, keyStreamByteDiffusion - Original keys
%   noiseType - 'gaussian' or 'salt & pepper'
%   noiseParam - Variance (gaussian) or Density (s&p)
%   dna_encode_func - Handle to @dna_encode
%   decryptImageDNA_func - Handle to @decryptImageDNA
% Output:
%   decryptedNoisyImage - The image after decrypting the noisy ciphertext
%   psnr_val - PSNR between originalImage and decryptedNoisyImage
%   mse_val - MSE between originalImage and decryptedNoisyImage

    fprintf('Running Noise Attack Analysis (%s)...\n', noiseType);
    originalSize = size(originalImage);
    numChannels = size(originalImage, 3);

    % --- 1. Add Noise to Encrypted Image ---
    try
        fprintf('  Adding %s noise (Param=%.4g) to encrypted image...\n', noiseType, noiseParam);
        noisyEncryptedUint8 = imnoise(encryptedImageUint8, noiseType, noiseParam);
    catch ME_imnoise
        error('Failed to add noise using imnoise: %s. Is Image Processing Toolbox available?', ME_imnoise.message);
    end
    % Display noisy encrypted image (optional)
    % figure; imshow(noisyEncryptedUint8); title(['Encrypted Image with ', noiseType, ' Noise']);

    % --- 2. Re-encode Noisy Uint8 back to DNA ---
    % This simulates the attacker having the noisy ciphertext, which needs
    % to be in the format expected by the decryption function (DNA sequence).
    fprintf('  Re-encoding noisy uint8 image to DNA sequence...\n');
    numPixelsPerChannel = numel(originalImage(:,:,1));
    try
        if numChannels == 1
            noisyEncryptedDna = dna_encode_func(noisyEncryptedUint8, keyStreamBinaryRule, numPixelsPerChannel);
        else % Color
            noisyEncryptedDna = cell(1, 3);
            for k = 1:3
                noisyEncryptedDna{k} = dna_encode_func(noisyEncryptedUint8(:,:,k), keyStreamBinaryRule, numPixelsPerChannel);
            end
        end
    catch ME_encode
         error('Failed to re-encode noisy image to DNA: %s', ME_encode.message);
    end
    fprintf('  Re-encoding complete.\n');

    % --- 3. Decrypt the Noisy DNA Sequence ---
    fprintf('  Decrypting the noisy DNA sequence...\n');
    try
        decryptedNoisyImage = decryptImageDNA_func(noisyEncryptedDna, originalSize, keyStreamBinaryRule, keyStreamByteDiffusion);
    catch ME_decrypt
         error('Failed to decrypt noisy DNA sequence: %s', ME_decrypt.message);
    end
     fprintf('  Decryption of noisy data complete.\n');

    % --- 4. Calculate PSNR/MSE ---
    fprintf('  Calculating PSNR/MSE against original image...\n');
    try
         if numChannels == 1
            mse_val = immse(decryptedNoisyImage, originalImage);
            psnr_val = psnr(decryptedNoisyImage, originalImage);
        else % Color: Average
            mse_vals = zeros(1,3); psnr_vals = zeros(1,3);
            for k=1:3
                mse_vals(k) = immse(decryptedNoisyImage(:,:,k), originalImage(:,:,k));
                 if mse_vals(k) < 1e-10
                     psnr_vals(k) = Inf; % Handle perfect channel recovery if noise was zero
                 else
                     psnr_vals(k) = psnr(decryptedNoisyImage(:,:,k), originalImage(:,:,k));
                 end
            end
            mse_val = mean(mse_vals);
            psnr_val = mean(psnr_vals(isfinite(psnr_vals))); % Average finite PSNRs
            if isempty(psnr_val), psnr_val = Inf; end % All channels perfect
        end
    catch ME_psnr
        warning('Could not calculate PSNR/MSE for noisy decryption: %s', ME_psnr.message);
        psnr_val = NaN;
        mse_val = NaN;
    end

    fprintf('Noise attack analysis finished.\n');
end
