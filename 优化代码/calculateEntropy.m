function entropy = calculateEntropy(imageChannel)
% Calculates the information entropy of a single image channel (grayscale).
% Input:
%   imageChannel: Single channel image (uint8 HxW)
% Output:
%   entropy: Calculated information entropy

    if ndims(imageChannel) ~= 2 || ~isa(imageChannel, 'uint8')
        error('Input must be a single-channel uint8 image (HxW).');
    end

    counts = imhist(imageChannel); % Get pixel value counts (0-255)
    totalPixels = numel(imageChannel);

    % Calculate probability P(i) = counts(i) / totalPixels
    probabilities = counts / totalPixels;

    % Remove zero probability entries to avoid log2(0)
    probabilities = probabilities(probabilities > 0);

    % Calculate entropy H = -sum(P(i) * log2(P(i)))
    entropy = -sum(probabilities .* log2(probabilities));
end

