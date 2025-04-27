function [npcr_avg, uaci_avg] = calculateNPCR_UACI(image1, image2)
% Calculates average NPCR and UACI between two images (Grayscale or Color).
% For color images, calculates per channel and averages the results.
% Input:
%   image1 - First image (uint8 HxW or HxWx3)
%   image2 - Second image (uint8, same size as image1)
% Output:
%   npcr_avg - Average NPCR across channels (%)
%   uaci_avg - Average UACI across channels (%)

    if ~isequal(size(image1), size(image2))
        error('Images must have the same dimensions for NPCR/UACI calculation.');
    end
    if ~isa(image1, 'uint8') || ~isa(image2, 'uint8')
         error('Input images must be uint8 type.');
    end

    imageSize = size(image1);
    rows = imageSize(1);
    cols = imageSize(2);
    numChannels = size(image1, 3);
    numPixelsPerChannel = rows * cols;

    npcr_vals = zeros(1, numChannels);
    uaci_vals = zeros(1, numChannels);

    for k = 1:numChannels
        channel1 = image1(:,:,k);
        channel2 = image2(:,:,k);

        % --- NPCR Calculation for channel k ---
        diff_pixels_matrix = (channel1 ~= channel2);
        diff_pixels_count = sum(diff_pixels_matrix(:));
        npcr_vals(k) = (diff_pixels_count / numPixelsPerChannel) * 100;

        % --- UACI Calculation for channel k ---
        diff_intensity = abs(double(channel1) - double(channel2));
        uaci_vals(k) = (sum(diff_intensity(:)) / (255 * numPixelsPerChannel)) * 100;
    end

    % --- Average Results ---
    npcr_avg = mean(npcr_vals);
    uaci_avg = mean(uaci_vals);

    % Optional: Print per-channel results if needed
    % if numChannels > 1
    %     fprintf('NPCR per channel: R=%.2f%%, G=%.2f%%, B=%.2f%%\n', npcr_vals(1), npcr_vals(2), npcr_vals(3));
    %     fprintf('UACI per channel: R=%.2f%%, G=%.2f%%, B=%.2f%%\n', uaci_vals(1), uaci_vals(2), uaci_vals(3));
    % end
end
