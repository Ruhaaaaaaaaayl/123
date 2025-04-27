function binaryDataString = dna_decode(dnaSequence, keyStreamBinaryRule, numPixels)
% Vectorized DNA decoding using key-stream selected rules.
% Input:
%   dnaSequence - DNA sequence string (char 1 x (numPixels * 4))
%   keyStreamBinaryRule - Binary key for rule selection (char 1 x (numPixels * 4))
%   numPixels - Number of pixels corresponding to the DNA sequence
% Output:
%   binaryDataString - Decoded binary data string (char 1 x (numPixels * 8))

    if ~(ischar(dnaSequence) || isstring(dnaSequence)) || length(dnaSequence) ~= numPixels * 4
        error('Input dnaSequence must be a char array of length numPixels*4.');
    end
    if ~(ischar(keyStreamBinaryRule) || isstring(keyStreamBinaryRule)) || length(keyStreamBinaryRule) ~= numPixels * 4
        error('Rule key stream must be a char array of length numPixels*4.');
    end

    % --- DNA Decoding Rules Lookup Table ---
    % Rows: Rules 0-7. Columns: Bases 'A', 'C', 'G', 'T' (indices 1-4).
    % Content: Binary pair value (0-3).
    decoding_table = [
        0 2 1 3; % Rule 0: A=00(0), C=10(2), G=01(1), T=11(3)
        0 3 1 2; % Rule 1: A=00(0), C=11(3), G=01(1), T=10(2)
        1 2 0 3; % Rule 2: A=01(1), C=10(2), G=00(0), T=11(3)
        1 3 0 2; % Rule 3: A=01(1), C=11(3), G=00(0), T=10(2)
        2 0 3 1; % Rule 4: A=10(2), C=00(0), G=11(3), T=01(1)
        3 0 2 1; % Rule 5: A=11(3), C=00(0), G=10(2), T=01(1)
        2 1 3 0; % Rule 6: A=10(2), C=01(1), G=11(3), T=00(0)
        3 1 2 0  % Rule 7: A=11(3), C=01(1), G=10(2), T=00(0)
    ]; % 8x4 double array

    % --- Prepare for Vectorized Lookup ---
    % 1. Map DNA bases 'A', 'C', 'G', 'T' to numeric column indices 1-4
    dnaChars = dnaSequence(:); % Ensure column vector (N*4) x 1 char
    dna_numeric_indices = zeros(size(dnaChars)); % (N*4) x 1 double
    dna_numeric_indices(dnaChars == 'A') = 1;
    dna_numeric_indices(dnaChars == 'C') = 2;
    dna_numeric_indices(dnaChars == 'G') = 3;
    dna_numeric_indices(dnaChars == 'T') = 4;
    if any(dna_numeric_indices == 0)
        error('DNA sequence contains invalid base characters.');
    end

    % 2. Determine rule index (0-7) for each pixel using first 3 bits of key
    keyStreamBinaryRuleMatrix = reshape(keyStreamBinaryRule, 4, numPixels)'; % N x 4 char
    rule_indices_binary = keyStreamBinaryRuleMatrix(:, 1:3); % N x 3 char
    rule_indices = bin2dec(rule_indices_binary); % N x 1 double (0-7)

    % 3. Expand rule indices for each base (row indices)
    rule_indices_expanded = repelem(rule_indices, 4); % (N*4) x 1 double

    % --- Vectorized Decoding using Linear Indexing ---
    % MATLAB indices are 1-based
    linear_indices = sub2ind(size(decoding_table), rule_indices_expanded + 1, dna_numeric_indices);
    binary_pair_values = decoding_table(linear_indices); % (N*4) x 1 double (0-3)

    % --- Convert binary pair values back to 2-bit binary strings ---
    binary_pairs_char = dec2bin(binary_pair_values, 2); % (N*4) x 2 char array

    % --- Reshape to final 8-bit binary string ---
    binaryDataString = reshape(binary_pairs_char', 1, numPixels * 8); % 1 x (N*8) char row vector
end
