% =========================================================================
% Main Script for Lorenz+DNA Color/Grayscale Image Encryption & Analysis
% =========================================================================
clear; clc; close all; % Clean workspace

disp('--- Lorenz+DNA ͼ���������� (��ɫ/�Ҷ�) ---'); % <-- (��΢�޸�������)

% --- ��������� ---
% ѡ��ͼ��
[fileName, pathName] = uigetfile({'*.png'; '*.jpg'; '*.jpeg'; '*.bmp'; '*.tif'; '*.tiff'}, 'ѡ��һ��ͼ���ļ� (��ɫ��Ҷ�)');
if isequal(fileName, 0) || isequal(pathName, 0)
    disp('�û�ȡ��ѡ�������˳���');
    return;
else
    imagePath = fullfile(pathName, fileName);
    disp(['ѡ���ͼ��: ', imagePath]);
end

% Lorenz ���� (����ֵ, ��������Կ��һ����)
rho = 28;
sigma = 10;
beta = 8/3;
% Lorenz ��ʼ���� (��Կ�Ĺؼ�����!)
initialConditions = [0.11, 0.22, 0.33]; % ʹ�÷�ƽ������ͬ��ֵ
key_sensitivity_delta = 1e-14; % ���������Բ��Ե�΢С�仯

% ��������
noise_type = 'gaussian'; % 'gaussian' (��˹), 'salt & pepper' (����)
noise_param_gaussian = 0.001; % ��˹�����ķ���
noise_param_sp = 0.05; % �����������ܶ�
crop_percentage = 0.1; % �����Ĳü��Ŀ��/�߶Ȱٷֱ� (����, 0.1 = 10%)

% --- 1. ��ȡԭʼͼ�񲢻�ȡ��Ϣ ---
try
    originalImage = imread(imagePath);
catch ME
    error('�޷���ȡͼ���ļ�: %s������·�����ļ������ԡ�\n%s', imagePath, ME.message);
end
originalImage = uint8(originalImage); % ȷ��Ϊ uint8 ����
originalSize = size(originalImage);
rows = originalSize(1);
cols = originalSize(2);
numChannels = size(originalImage, 3); % �Ҷ�ͼΪ1, ��ɫͼΪ3
numPixelsPerChannel = rows * cols;
totalPixels = numel(originalImage);

fprintf('ͼ��ߴ�: %d x %d x %d (%d ������)\n', rows, cols, numChannels, totalPixels);
if numChannels == 3
    disp('�����ǲ�ɫͼ�� (RGB)��');
elseif numChannels == 1
    disp('�����ǻҶ�ͼ��');
else
    error('��֧�ֵ�ͼ���ʽ��ͨ�������� 1 �� 3��');
end

% --- 2. ���� Lorenz ��Կ�� ---
% ��Կ����ÿͨ����С���ɣ���Ϊ R, G, B ����
fprintf('Ϊ %d �������� Lorenz ��Կ��...\n', numPixelsPerChannel);
tic;
try
    [keyStreamBinaryRule, keyStreamByteDiffusion] = generateLorenzKeyStream(initialConditions, rho, sigma, beta, numPixelsPerChannel);
catch ME
    fprintf('��Կ���ɹ����г���: %s\n', ME.message);
    fprintf('�ļ�: %s, ��: %d\n', ME.stack(1).file, ME.stack(1).line);
    return;
end
keyGenTime = toc;
fprintf('��Կ�����ɳɹ� (%.4f ��)��\n', keyGenTime);

% --- 3. ����ͼ�� ---
fprintf('����ͼ�� (%d ͨ��)...\n', numChannels);
tic;
try
    % ��������ͼ�� (��ɫ��Ҷ�)
    [encryptedDnaData, sizeCheck] = encryptImageDNA(originalImage, keyStreamBinaryRule, keyStreamByteDiffusion);
    % encryptedDnaData ���ڲ�ɫ�� cell ����, ���ڻҶ��� char ����
catch ME
    fprintf('���ܹ����г���: %s\n', ME.message);
    fprintf('�ļ�: %s, ��: %d\n', ME.stack(1).file, ME.stack(1).line);
    return;
end
encryptTime = toc;
fprintf('������� (%.4f ��)��\n', encryptTime);
if ~isequal(sizeCheck, originalSize)
    warning('���ܺ������صĳߴ粻ƥ�䣡');
end

% --- 4. ��ȡ���ڷ���/��ʾ�� uint8 ����ͼ�� ---
% ���� DNA (����������ɢ) �Ի�� uint8 ��ʾ
fprintf('�������ڷ����� uint8 ����ͼ��...\n');
encryptedImageUint8 = zeros(originalSize, 'uint8');
try
    if numChannels == 1
        decodedBinaryString = dna_decode(encryptedDnaData, keyStreamBinaryRule, numPixelsPerChannel);
        decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
        encryptedImageUint8(:) = uint8(bin2dec(decodedBinaryMatrix));
        encryptedImageUint8 = reshape(encryptedImageUint8, rows, cols);
    else % Color
        for k = 1:numChannels
            decodedBinaryString = dna_decode(encryptedDnaData{k}, keyStreamBinaryRule, numPixelsPerChannel);
            decodedBinaryMatrix = reshape(decodedBinaryString, 8, numPixelsPerChannel)';
            channel_uint8 = uint8(bin2dec(decodedBinaryMatrix));
            encryptedImageUint8(:,:,k) = reshape(channel_uint8, rows, cols);
        end
    end
catch ME
    fprintf('���� uint8 ����ͼ��ʱ����: %s\n', ME.message);
    return;
end
fprintf('Uint8 ����ͼ�������ɡ�\n');

% --- 5. ��ʾͼ�� (ԭʼ, ����) ---
figure('Name', '���ܽ��', 'NumberTitle', 'off'); % <--- �޸�����
subplot(1, 2, 1); imshow(originalImage); title('ԭʼͼ��'); % <--- �޸�����
subplot(1, 2, 2); imshow(encryptedImageUint8); title('����ͼ�� (uint8)'); % <--- �޸�����

% --- 6. ֱ��ͼ���� ---
figure('Name', 'ֱ��ͼ����', 'NumberTitle', 'off'); % <-- (��ѡ) �޸� Figure Name
if numChannels == 1
    subplot(1, 2, 1); imhist(originalImage); title('ԭʼֱ��ͼ'); ylim('auto'); % <-- (��ѡ) �޸� title
    subplot(1, 2, 2); imhist(encryptedImageUint8); title('����ֱ��ͼ'); ylim('auto'); % <-- (��ѡ) �޸� title
else % Color: Show histograms for each channel
    colorLabels = {'��ɫ', '��ɫ', '��ɫ'}; % <-- (��ѡ) �޸ı�ǩ
    for k = 1:3
        subplot(2, 3, k);
        imhist(originalImage(:,:,k)); title(['ԭʼ', colorLabels{k}, 'ͨ��']); ylim('auto'); % <-- (��ѡ) �޸� title
        subplot(2, 3, k+3);
        imhist(encryptedImageUint8(:,:,k)); title(['����', colorLabels{k}, 'ͨ��']); ylim('auto'); % <-- (��ѡ) �޸� title
    end
    % Adjust layout slightly
    annotation('textbox', [0 0.9 1 0.1], 'String', 'ֱ��ͼ���� (RGB ͨ��)', ... % <-- (��ѡ) �޸� Annotation
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% --- 7. ��Ϣ�ط��� ---
fprintf('\n--- ��Ϣ�ط��� ---\n');
totalEntropyOrig = 0;
totalEntropyEnc = 0;
if numChannels == 1
    entropyOrig = calculateEntropy(originalImage);
    entropyEnc = calculateEntropy(encryptedImageUint8);
    fprintf('ԭʼ�Ҷ���: %.4f\n', entropyOrig);
    fprintf('���ܻҶ���: %.4f (����ֵ ~8)\n', entropyEnc);
    totalEntropyOrig = entropyOrig;
    totalEntropyEnc = entropyEnc;
else
    fprintf('ͨ����:\n');
    fprintf('%10s %10s %10s\n', 'ͨ��', 'ԭʼ', '����');
    colorLabels = {'��', '��', '��'}; % ʹ�����涨������ı�ǩ
    for k = 1:3
        entropyOrig = calculateEntropy(originalImage(:,:,k));
        entropyEnc = calculateEntropy(encryptedImageUint8(:,:,k));
        fprintf('%10s %10.4f %10.4f\n', colorLabels{k}, entropyOrig, entropyEnc);
        totalEntropyOrig = totalEntropyOrig + entropyOrig;
        totalEntropyEnc = totalEntropyEnc + entropyEnc;
    end
    fprintf('ƽ��ԭʼ��: %.4f\n', totalEntropyOrig / 3);
    fprintf('ƽ��������: %.4f (����ֵ ~8)\n', totalEntropyEnc / 3);
end

% --- 8. ����Է��� ---
fprintf('\n--- ������������Է��� ---\n');
% �Բ�ɫͼ���ÿ��ͨ�����з���
analyzeCorrelation(originalImage, 'ԭʼͼ��'); % <-- (��ѡ) �޸ı�ǩ
analyzeCorrelation(encryptedImageUint8, '����ͼ��'); % ���� uint8 ���ܰ汾, <-- (��ѡ) �޸ı�ǩ

% --- 9. NPCR & UACI ���� ---
fprintf('\n--- NPCR & UACI ���� ---\n');
try
    % ����һ����΢�޸ĵ�ԭʼͼ�� (��תһ�������е�һ������)
    originalImage_modified = originalImage;
    rand_row = randi([1, rows]);
    rand_col = randi([1, cols]);
    rand_channel = randi([1, numChannels]); % ѡ��һ��ͨ�������޸�
    originalImage_modified(rand_row, rand_col, rand_channel) = bitxor(originalImage(rand_row, rand_col, rand_channel), 1);

    % ʹ����ͬ����Կ�����޸ĺ��ͼ��
    fprintf('Ϊ NPCR/UACI �����޸ĺ��ͼ��...\n');
    [encryptedDnaData_modified, ~] = encryptImageDNA(originalImage_modified, keyStreamBinaryRule, keyStreamByteDiffusion);

    % �����޸ĺ����ͼ��� uint8 �汾
    encryptedImageUint8_modified = zeros(originalSize, 'uint8');
     if numChannels == 1
        decodedBinaryString_mod = dna_decode(encryptedDnaData_modified, keyStreamBinaryRule, numPixelsPerChannel);
        decodedBinaryMatrix_mod = reshape(decodedBinaryString_mod, 8, numPixelsPerChannel)';
        encryptedImageUint8_modified(:) = uint8(bin2dec(decodedBinaryMatrix_mod));
        encryptedImageUint8_modified = reshape(encryptedImageUint8_modified, rows, cols);
    else % Color
        for k = 1:numChannels
            decodedBinaryString_mod = dna_decode(encryptedDnaData_modified{k}, keyStreamBinaryRule, numPixelsPerChannel);
            decodedBinaryMatrix_mod = reshape(decodedBinaryString_mod, 8, numPixelsPerChannel)';
            channel_uint8_mod = uint8(bin2dec(decodedBinaryMatrix_mod));
            encryptedImageUint8_modified(:,:,k) = reshape(channel_uint8_mod, rows, cols);
        end
     end
    fprintf('�޸ĺ�ļ���ͼ�������ɡ�\n');

    % ���� NPCR �� UACI (�����ڲ�������ɫ)
    [npcr_avg, uaci_avg] = calculateNPCR_UACI(encryptedImageUint8, encryptedImageUint8_modified);
    fprintf('ƽ�� NPCR: %.4f%% (����ֵ > 99.6%%)\n', npcr_avg);
    fprintf('ƽ�� UACI: %.4f%% (���� uint8 ����ֵ ~ 33.4%%)\n', uaci_avg);

catch ME
    fprintf('NPCR/UACI ��������г���: %s\n', ME.message);
    fprintf('�ļ�: %s, ��: %d\n', ME.stack(1).file, ME.stack(1).line);
end

% --- 10. ���� ---
fprintf('\n����ͼ��...\n');
tic;
try
    % ���ݼ��ܵ� DNA ���� (cell �� char) ��ԭʼ�ߴ�
    decryptedImage = decryptImageDNA(encryptedDnaData, originalSize, keyStreamBinaryRule, keyStreamByteDiffusion);
catch ME
    fprintf('���ܹ����г���: %s\n', ME.message);
    fprintf('�ļ�: %s, ��: %d\n', ME.stack(1).file, ME.stack(1).line);
    decryptedImage = []; % ��ǽ���ʧ��
end
decryptTime = toc;

if ~isempty(decryptedImage)
    fprintf('������� (%.4f ��)��\n', decryptTime);

    % --- 11. ��ʾ���ܽ������֤ ---
    figure('Name', '���ܽ������֤', 'NumberTitle', 'off'); % <--- �޸�����
    subplot(1, 2, 1); imshow(originalImage); title('ԭʼͼ��'); % <--- �޸�����
    subplot(1, 2, 2); imshow(decryptedImage); title('����ͼ��'); % <--- �޸�����

    % ���� MSE �� PSNR
    fprintf('\n--- ������������ ---\n');
    try
        if numChannels == 1
            mse = immse(decryptedImage, originalImage);
            psnr_val = psnr(decryptedImage, originalImage);
        else % Color: �����ͨ����ƽ�� MSE/PSNR
            mse_r = immse(decryptedImage(:,:,1), originalImage(:,:,1));
            mse_g = immse(decryptedImage(:,:,2), originalImage(:,:,2));
            mse_b = immse(decryptedImage(:,:,3), originalImage(:,:,3));
            mse = mean([mse_r, mse_g, mse_b]);

            % PSNR ��Ҫ���� MSE�����������ؽ��������
            if mse < 1e-10
                 psnr_val = Inf;
            else
                psnr_r = psnr(decryptedImage(:,:,1), originalImage(:,:,1));
                psnr_g = psnr(decryptedImage(:,:,2), originalImage(:,:,2));
                psnr_b = psnr(decryptedImage(:,:,3), originalImage(:,:,3));
                % �� PSNR dB ֵ��ƽ���ǳ�������
                psnr_val = mean([psnr_r, psnr_g, psnr_b]);
            end
        end

        fprintf('������� (MSE): %.4g\n', mse);
        if isinf(psnr_val) || psnr_val > 100 % ��� Inf ��ǳ��ߵ� PSNR
             fprintf('��ֵ����� (PSNR): Inf dB\n');
             fprintf('�ɹ�: ����ͼ����ԭʼͼ����ȫƥ�䣡\n'); % <-- (��ѡ) �޸��ı�
        else
            fprintf('��ֵ����� (PSNR): %.4f dB\n', psnr_val);
            if mse > 1e-6 % ����΢С�ĸ������
                 fprintf('����: ����ͼ����ԭʼͼ�����в�ͬ��\n'); % <-- (��ѡ) �޸��ı�
            else
                 fprintf('�ɹ�: ����ͼ�����Ӿ�����ԭʼͼ����ͬ��\n'); % <-- (��ѡ) �޸��ı�
            end
        end
    catch ME_psnr
        fprintf('�޷����� PSNR/MSE: %s\n', ME_psnr.message);
    end
else
    fprintf('����ʧ�ܡ��޷�����������\n');
end

% --- 12. ��Կ�ռ���� ---
fprintf('\n--- ��Կ�ռ���� ---\n');
analyzeKeySpace(initialConditions, [rho, sigma, beta]);

% --- 13. ��Կ�����Է��� ---
fprintf('\n--- ��Կ�����Է��� ---\n');
try
    % �����������ĺ��ĺ������
    encrypt_func_handle = @encryptImageDNA;
    decrypt_func_handle = @decryptImageDNA; % ���ﲻ���ϸ���Ҫ��������ʵ��
    keygen_func_handle = @generateLorenzKeyStream;
    dna_decode_func_handle = @dna_decode; % ��Ҫ��ȡ uint8 ���ڱȽ�
    npcr_uaci_func_handle = @calculateNPCR_UACI;

    [npcr_sens, uaci_sens] = analyzeKeySensitivity(...
        originalImage, ...
        encrypt_func_handle, ...
        keygen_func_handle, ...
        dna_decode_func_handle, ...
        npcr_uaci_func_handle, ...
        initialConditions, rho, sigma, beta, ...
        key_sensitivity_delta);

    fprintf('����Կ�仯�������� (delta=%.0e):\n', key_sensitivity_delta);
    fprintf('  ����֮��� NPCR: %.4f%%\n', npcr_sens);
    fprintf('  ����֮��� UACI: %.4f%%\n', uaci_sens);
    if npcr_sens > 99 && uaci_sens > 30
        fprintf('  ���: ��⵽����Կ������ (���ã�)��\n'); % <-- (��ѡ) �޸��ı�
    else
        fprintf('  ���: ��⵽����Կ������ (Ǳ�����㣡)��\n'); % <-- (��ѡ) �޸��ı�
    end
catch ME_sens
    fprintf('��Կ�����Է��������г���: %s\n', ME_sens.message);
    fprintf('�ļ�: %s, ��: %d\n', ME_sens.stack(1).file, ME_sens.stack(1).line);
end


% --- 14. ������������ ---
fprintf('\n--- ������������ ---\n');
if ~isempty(decryptedImage) % ������ʼ����/���ܳɹ�ʱ����
    try
        % ��������ѡ����������
        noise_type_cn = ''; % ���ڱ������������
        if strcmpi(noise_type, 'gaussian')
            noise_param = noise_param_gaussian;
            noise_type_cn = '��˹';
            fprintf('Ӧ�ø�˹���� (���� = %.4f)...\n', noise_param);
        elseif strcmpi(noise_type, 'salt & pepper')
            noise_param = noise_param_sp;
             noise_type_cn = '����';
            fprintf('Ӧ�ý������� (�ܶ� = %.2f)...\n', noise_param);
        else
            error('��֧�ֵ���������: %s', noise_type);
        end

        % ����Ҫ������ݸ���������
        [decryptedNoisyImage, psnr_noise, mse_noise] = analyzeNoiseAttack(...
            originalImage, ...
            encryptedImageUint8, ... % ���� uint8 �汾
            keyStreamBinaryRule, ...
            keyStreamByteDiffusion, ...
            noise_type, noise_param, ...
            @dna_encode, ... % ���ݺ������
            @decryptImageDNA); % ���ݺ������

        fprintf('��������������ɡ�\n');
        fprintf('  MSE (ԭʼ vs. ���ܵĺ���ͼ��): %.4f\n', mse_noise);
        fprintf('  PSNR (ԭʼ vs. ���ܵĺ���ͼ��): %.4f dB\n', psnr_noise);

        figure('Name', '�����������', 'NumberTitle', 'off'); % <-- (��ѡ) �޸� Figure Name
        subplot(1, 2, 1); imshow(originalImage); title('ԭʼͼ��'); % <-- ʹ�����ı���
        subplot(1, 2, 2); imshow(decryptedNoisyImage); title(sprintf('%s�������������', noise_type_cn)); % <--- �޸�����

    catch ME_noise
        fprintf('�����������������г���: %s\n', ME_noise.message);
         fprintf('�ļ�: %s, ��: %d\n', ME_noise.stack(1).file, ME_noise.stack(1).line);
    end
else
    fprintf('���ڳ�ʼ����ʧ�ܣ�������������������\n');
end


% --- 15. �ü��������� ---
fprintf('\n--- �ü��������� ---\n');
if ~isempty(decryptedImage) % ������ʼ����/���ܳɹ�ʱ����
    try
        fprintf('Ӧ�òü����� (�����Ĳü� %.1f%%)...\n', crop_percentage*100);

        % ����ü����� (����)
        crop_width = floor(cols * crop_percentage);
        crop_height = floor(rows * crop_percentage);
        start_x = floor((cols - crop_width) / 2) + 1;
        start_y = floor((rows - crop_height) / 2) + 1;
        % ������������ͼ��߽���
        start_x = max(1, start_x);
        start_y = max(1, start_y);
        end_x = min(cols, start_x + crop_width - 1);
        end_y = min(rows, start_y + crop_height - 1);

        cropRectangle = [start_x, start_y, end_x - start_x + 1, end_y - start_y + 1]; % ĳЩ���������� [x, y, width, height] ��ʽ

        % ����Ҫ������ݸ���������
        [decryptedCroppedImage, psnr_crop, mse_crop] = analyzeCroppingAttack(...
            originalImage, ...
            encryptedImageUint8, ... % ���� uint8 �汾
            keyStreamBinaryRule, ...
            keyStreamByteDiffusion, ...
            cropRectangle, ... % ���ݼ�����ľ���
            @dna_encode, ... % ���ݺ������
            @decryptImageDNA); % ���ݺ������

         fprintf('�ü�����������ɡ�\n');
         fprintf('  MSE (ԭʼ vs. ���ܵĲü�ͼ��): %.4f\n', mse_crop);
         fprintf('  PSNR (ԭʼ vs. ���ܵĲü�ͼ��): %.4f dB\n', psnr_crop);

        figure('Name', '�ü��������', 'NumberTitle', 'off'); % <-- (��ѡ) �޸� Figure Name
        subplot(1, 2, 1); imshow(originalImage); title('ԭʼͼ��'); % <-- ʹ�����ı���
        subplot(1, 2, 2); imshow(decryptedCroppedImage); title(sprintf('%.0f%% �ü����������', crop_percentage*100)); % <--- �޸�����

    catch ME_crop
        fprintf('�ü��������������г���: %s\n', ME_crop.message);
        fprintf('�ļ�: %s, ��: %d\n', ME_crop.stack(1).file, ME_crop.stack(1).line);
    end
else
     fprintf('���ڳ�ʼ����ʧ�ܣ������ü�����������\n');
end

fprintf('\n--- ���д���ͷ������ ---\n');
