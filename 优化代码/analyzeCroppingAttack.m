function [decryptedCroppedImage, psnr_val, mse_val] = analyzeCroppingAttack(...
            originalImage, encryptedImageUint8, ...
            keyStreamBinaryRule, keyStreamByteDiffusion, ...
            cropRectangle, ... % [x_start, y_start, width, height]
            dna_encode_func, decryptImageDNA_func)
% Simulates a cropping attack on the encrypted image and decrypts it.
% Input:
%   originalImage - The original clean image (for PSNR comparison)
%   encryptedImageUint8 - The encrypted image (uint8 HxW or HxWx3) to attack
%   keyStreamBinaryRule, keyStreamByteDiffusion - Original keys
%   cropRectangle - [x_start, y_start, width, height] defining the crop area
%   dna_encode_func - Handle to @dna_encode
%   decryptImageDNA_func - Handle to @decryptImageDNA
% Output:
%   decryptedCroppedImage - Image after decrypting the cropped ciphertext
%   psnr_val - PSNR between originalImage and decryptedCroppedImage
%   mse_val - MSE between originalImage and decryptedCroppedImage

    fprintf('Running Cropping Attack Analysis...\n');
    originalSize = size(originalImage);
    numChannels = size(originalImage, 3);
    rows = originalSize(1);
    cols = originalSize(2);

    % --- 1. Apply Cropping to Encrypted Image ---
    % Simulate cropping by setting the cropped area to a constant value (e.g., 0 or 128)
    % This represents data loss in the ciphertext.
    croppedEncryptedUint8 = encryptedImageUint8;
    x1 = floor(cropRectangle(1));
    y1 = floor(cropRectangle(2));
    w = floor(cropRectangle(3));
    h = floor(cropRectangle(4));
    x2 = min(cols, x1 + w - 1); % Ensure coordinates are within bounds
    y2 = min(rows, y1 + h - 1);
    x1 = max(1, x1);
    y1 = max(1, y1);

    cropValue = 0; % Value to fill the cropped area (e.g., 0 for black)
    fprintf('  Applying cropping (setting region [%d:%d, %d:%d] to %d) to encrypted image...\n', y1, y2, x1, x2, cropValue);

    if y1 <= y2 && x1 <= x2 % Check if rectangle has valid dimensions
      croppedEncryptedUint8(y1:y2, x1:x2, :) = cropValue;
    else
      warning('Cropping rectangle has zero or negative size, skipping crop application.');
    end

    % Display cropped encrypted image (optional)
    % figure; imshow(croppedEncryptedUint8); title('Encrypted Image after Cropping');

    % --- 2. Re-encode Cropped Uint8 back to DNA ---
    fprintf('  Re-encoding cropped uint8 image to DNA sequence...\n');
    numPixelsPerChannel = rows * cols;
     try
        if numChannels == 1
            croppedEncryptedDna = dna_encode_func(croppedEncryptedUint8, keyStreamBinaryRule, numPixelsPerChannel);
        else % Color
            croppedEncryptedDna = cell(1, 3);
            for k = 1:3
                croppedEncryptedDna{k} = dna_encode_func(croppedEncryptedUint8(:,:,k), keyStreamBinaryRule, numPixelsPerChannel);
            end
        end
     catch ME_encode
         error('Failed to re-encode cropped image to DNA: %s', ME_encode.message);
     end
    fprintf('  Re-encoding complete.\n');

    % --- 3. Decrypt the Cropped DNA Sequence ---
    fprintf('  Decrypting the cropped DNA sequence...\n');
    try
        decryptedCroppedImage = decryptImageDNA_func(croppedEncryptedDna, originalSize, keyStreamBinaryRule, keyStreamByteDiffusion);
    catch ME_decrypt
         error('Failed to decrypt cropped DNA sequence: %s', ME_decrypt.message);
    end
     fprintf('  Decryption of cropped data complete.\n');

    % --- 4. Calculate PSNR/MSE ---
    fprintf('  Calculating PSNR/MSE against original image...\n');
    try
         if numChannels == 1
            mse_val = immse(decryptedCroppedImage, originalImage);
            psnr_val = psnr(decryptedCroppedImage, originalImage);
        else % Color: Average
            mse_vals = zeros(1,3); psnr_vals = zeros(1,3);
             for k=1:3
                mse_vals(k) = immse(decryptedCroppedImage(:,:,k), originalImage(:,:,k));
                 if mse_vals(k) < 1e-10
                     psnr_vals(k) = Inf;
                 else
                     psnr_vals(k) = psnr(decryptedCroppedImage(:,:,k), originalImage(:,:,k));
                 end
            end
            mse_val = mean(mse_vals);
            psnr_val = mean(psnr_vals(isfinite(psnr_vals)));
            if isempty(psnr_val), psnr_val = Inf; end
        end
    catch ME_psnr
        warning('Could not calculate PSNR/MSE for cropped decryption: %s', ME_psnr.message);
        psnr_val = NaN;
        mse_val = NaN;
    end

    fprintf('Cropping attack analysis finished.\n');
end
