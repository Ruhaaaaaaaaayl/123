function dnaSequence = dna_encode(inputData, keyStreamBinaryRule, numPixels)
% Vectorized DNA encoding using key-stream selected rules.
% Input:
%   inputData - uint8 image vector (1 x numPixels) or matrix (auto-flattened)
%   keyStreamBinaryRule - Binary key for rule selection (char 1 x (numPixels * 4))
%   numPixels - Number of pixels (length of inputData)
% Output:
%   dnaSequence - DNA sequence string (char 1 x (numPixels * 4))

    if ~isa(inputData, 'uint8')
        error('Input data must be uint8 type.');
    end
    if numel(inputData) ~= numPixels
       error('Input data size (%d) does not match numPixels (%d).', numel(inputData), numPixels);
    end
    if ~(ischar(keyStreamBinaryRule) || isstring(keyStreamBinaryRule)) || length(keyStreamBinaryRule) ~= numPixels * 4
        error('Rule key stream must be a char array of length numPixels*4.');
    end

    % --- Ensure row vector and convert to binary ---
    img_vector = inputData(:)'; % Flatten to 1 x N uint8 row vector
    binaryDataCharMatrix = dec2bin(img_vector, 8); % N x 8 char array

    % --- DNA Encoding Rules Lookup Table ---
    encoding_rules = [
        'A','G','C','T'; % Rule 0
        'A','G','T','C'; % Rule 1
        'G','A','C','T'; % Rule 2
        'G','A','T','C'; % Rule 3
        'C','T','A','G'; % Rule 4
        'C','T','G','A'; % Rule 5
        'T','C','A','G'; % Rule 6
        'T','C','G','A'  % Rule 7
    ]; % 8x4 char array

    % --- Prepare for Vectorized Lookup ---
    % 1. Reshape N x 8 binary matrix to (N*4) x 2 pairs
    binaryDataPairs = reshape(binaryDataCharMatrix', 2, numPixels * 4)'; % (N*4) x 2 char

    % 2. Convert '00'-'11' pairs to numeric values 0-3 (column indices)
    binary_pair_values = bin2dec(binaryDataPairs); % (N*4) x 1 double

    % 3. Determine rule index (0-7) for each pixel using first 3 bits of key
    keyStreamBinaryRuleMatrix = reshape(keyStreamBinaryRule, 4, numPixels)'; % N x 4 char
    rule_indices_binary = keyStreamBinaryRuleMatrix(:, 1:3); % N x 3 char
    rule_indices = bin2dec(rule_indices_binary); % N x 1 double (0-7)

    % 4. Expand rule indices for each base pair (row indices)
    rule_indices_expanded = repelem(rule_indices, 4); % (N*4) x 1 double

    % --- Vectorized Encoding using Linear Indexing ---
    % MATLAB indices are 1-based
    linear_indices = sub2ind(size(encoding_rules), rule_indices_expanded + 1, binary_pair_values + 1);
    dnaSequenceChars = encoding_rules(linear_indices); % (N*4) x 1 char

    % --- Format Output ---
    dnaSequence = dnaSequenceChars'; % Transpose to 1 x (N*4) char row vector
end
