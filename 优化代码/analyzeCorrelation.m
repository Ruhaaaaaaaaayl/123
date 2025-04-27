% =========================================================================
% Function for Pixel Correlation Analysis
% =========================================================================
function [corr_h, corr_v, corr_d] = analyzeCorrelation(imageData, imageType, colorChannel, numSamples)
% Performs correlation analysis on adjacent pixels (horizontal, vertical, diagonal).
%
% Args:
%   imageData: The image matrix (grayscale or single color channel).
%   imageType: String ('Original' or 'Encrypted') for plot titles.
%   colorChannel: String ('Grayscale', 'Red', 'Green', 'Blue') for plot titles.
%   numSamples: Number of pixel pairs to sample for correlation.
%
% Returns:
%   corr_h: Horizontal correlation coefficient.
%   corr_v: Vertical correlation coefficient.
%   corr_d: Diagonal correlation coefficient.
%
% Example Usage:
%   [ch, cv, cd] = CorrelationAnalysis(originalImage(:,:,1), 'Original', 'Red', 5000);

if nargin < 4
    numSamples = 3000; % Default number of samples if not provided
end
if nargin < 3
    colorChannel = 'Grayscale'; % Default channel if not provided
end
 if nargin < 2
     imageType = 'ͼ��'; % Default type if not provided
 end

fprintf('--- Correlation Analysis (%s - %s Channel) ---\n', imageType, colorChannel);

[height, width] = size(imageData);

% --- Translate Titles ---
imageType_CN = imageType; % Keep original input for potential logic checks
if strcmpi(imageType, 'Original')
    imageType_CN = 'ԭʼͼ��';
elseif strcmpi(imageType, 'Encrypted')
    imageType_CN = '����ͼ��';
end

colorChannel_CN = colorChannel; % Keep original input
switch lower(colorChannel)
    case 'red'
        colorChannel_CN = '��ɫ';
    case 'green'
        colorChannel_CN = '��ɫ';
    case 'blue'
        colorChannel_CN = '��ɫ';
    case 'grayscale'
        colorChannel_CN = '�Ҷ�';
    % Add other cases if needed
end


% --- Randomly Select Pixel Pairs ---
% Generate random indices, ensuring they are within bounds for adjacent pixels
rand_x = randi([1 height-1], numSamples, 1);
rand_y = randi([1 width-1], numSamples, 1);

% --- Extract Pixel Values ---
% Horizontal pairs: (x, y) and (x, y+1)
adj_h = zeros(numSamples, 2);
for k = 1:numSamples
    adj_h(k, 1) = imageData(rand_x(k), rand_y(k));
    adj_h(k, 2) = imageData(rand_x(k), rand_y(k) + 1);
end

% Vertical pairs: (x, y) and (x+1, y)
adj_v = zeros(numSamples, 2);
for k = 1:numSamples
    adj_v(k, 1) = imageData(rand_x(k), rand_y(k));
    adj_v(k, 2) = imageData(rand_x(k) + 1, rand_y(k));
end

% Diagonal pairs: (x, y) and (x+1, y+1)
adj_d = zeros(numSamples, 2);
for k = 1:numSamples
    adj_d(k, 1) = imageData(rand_x(k), rand_y(k));
    adj_d(k, 2) = imageData(rand_x(k) + 1, rand_y(k) + 1);
end

% --- Calculate Correlation Coefficients ---
corr_h = corr2(adj_h(:,1), adj_h(:,2));
corr_v = corr2(adj_v(:,1), adj_v(:,2));
corr_d = corr2(adj_d(:,1), adj_d(:,2));

fprintf('  Horizontal Correlation: %.4f\n', corr_h);
fprintf('  Vertical Correlation:   %.4f\n', corr_v);
fprintf('  Diagonal Correlation:   %.4f\n', corr_d);

% --- Plotting ---
figure('Name', sprintf('%s ����Է��� (%sͨ��)', imageType_CN, colorChannel_CN));

% Horizontal Plot
subplot(2, 2, 1);
plot(adj_h(:,1), adj_h(:,2), '.b');
xlabel('����ֵ p(i, j)');
ylabel('����ֵ p(i, j+1)');
% Combine title lines using cell array for multi-line title
title({sprintf('%s - %s (%d ����)', imageType_CN, colorChannel_CN, numSamples), ...
       sprintf('ˮƽ���� (���ϵ�� = %.4f)', corr_h)}); % <<< MODIFIED HERE
axis([0 255 0 255]); grid on; axis square;

% Vertical Plot
subplot(2, 2, 2);
plot(adj_v(:,1), adj_v(:,2), '.r');
xlabel('����ֵ p(i, j)');
ylabel('����ֵ p(i+1, j)');
title(sprintf('��ֱ���� (���ϵ�� = %.4f)', corr_v)); % <<< MODIFIED HERE
axis([0 255 0 255]); grid on; axis square;

% Diagonal Plot
subplot(2, 2, 3);
plot(adj_d(:,1), adj_d(:,2), '.g');
xlabel('����ֵ p(i, j)');
ylabel('����ֵ p(i+1, j+1)');
title(sprintf('�Խ��߷��� (���ϵ�� = %.4f)', corr_d)); % <<< MODIFIED HERE
axis([0 255 0 255]); grid on; axis square;

% Optional: Add a main title to the figure
% sgtitle(sprintf('%sͼ�������Է��� (%s ͨ��)', imageType_CN, colorChannel_CN)); % Requires R2018b or later

fprintf('--- Correlation plots generated for %s - %s ---\n', imageType, colorChannel);

end

