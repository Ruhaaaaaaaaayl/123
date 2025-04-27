% =========================================================================
% Main Script for Lorenz+DNA Color/Grayscale Image Encryption & Analysis
% =========================================================================
clear; clc; close all; % Clean workspace

disp('--- Lorenz+DNA 图像加密与分析 (彩色/灰度) ---'); % <-- (稍微修改了这里)

% --- 配置与参数 ---
% 选择图像
[fileName, pathName] = uigetfile({'*.png'; '*.jpg'; '*.jpeg'; '*.bmp'; '*.tif'; '*.tiff'}, '选择一个图像文件 (彩色或灰度)');
if isequal(fileName, 0) || isequal(pathName, 0)
    disp('用户取消选择。正在退出。');
    return;
else
    imagePath = fullfile(pathName, fileName);
    disp(['选择的图像: ', imagePath]);
end

% Lorenz 参数 (经典值, 可以是密钥的一部分)
rho = 28;
sigma = 10;
beta = 8/3;
% Lorenz 初始条件 (密钥的关键部分!)
initialConditions = [0.11, 0.22, 0.33]; % 使用非平凡、不同的值
key_sensitivity_delta = 1e-14; % 用于敏感性测试的微小变化

% 攻击参数
noise_type = 'gaussian'; % 'gaussian' (高斯), 'salt & pepper' (椒盐)
noise_param_gaussian = 0.001; % 高斯噪声的方差
noise_param_sp = 0.05; % 椒盐噪声的密度
crop_percentage = 0.1; % 从中心裁剪的宽度/高度百分比 (例如, 0.1 = 10%)

% --- 1. 读取原始图像并获取信息 ---
try
    originalImage = imread(imagePath);
catch ME
    error('无法读取图像文件: %s。请检查路径和文件完整性。\n%s', imagePath, ME.message);
end
originalImage = uint8(originalImage); % 确保为 uint8 类型
originalSize = size(originalImage);
rows = originalSize(1);
cols = originalSize(2);
numChannels = size(originalImage, 3); % 灰度图为1, 彩色图为3
numPixelsPerChannel = rows * cols;
totalPixels = numel(originalImage);

fprintf('图像尺寸: %d x %d x %d (%d 总像素)\n', rows, cols, numChannels, totalPixels);
if numChannels == 3
    disp('输入是彩色图像 (RGB)。');
elseif numChannels == 1
    disp('输入是灰度图像。');
else
    error('不支持的图像格式：通道数不是 1 或 3。');
end

% --- 2. 生成 Lorenz 密钥流 ---
% 密钥基于每通道大小生成，并为 R, G, B 复用
fprintf('为 %d 像素生成 Lorenz 密钥流...\n', numPixelsPerChannel);
tic;
try
    [keyStreamBinaryRule, keyStreamByteDiffusion] = generateLorenzKeyStream(initialConditions, rho, sigma, beta, numPixelsPerChannel);
catch ME
    fprintf('密钥生成过程中出错: %s\n', ME.message);
    fprintf('文件: %s, 行: %d\n', ME.stack(1).file, ME.stack(1).line);
    return;
end
keyGenTime = toc;
fprintf('密钥流生成成功 (%.4f 秒)。\n', keyGenTime);

% --- 3. 加密图像 ---
fprintf('加密图像 (%d 通道)...\n', numChannels);
tic;
try
    % 传递完整图像 (彩色或灰度)
    [encryptedDnaData, sizeCheck] = encryptImageDNA(originalImage, keyStreamBinaryRule, keyStreamByteDiffusion);
    % encryptedDnaData 对于彩色是 cell 数组, 对于灰度是 char 数组
catch ME
    fprintf('加密过程中出错: %s\n', ME.message);
    fprintf('文件: %s, 行: %d\n', ME.stack(1).file, ME.stack(1).line);
    return;
end
encryptTime = toc;
fprintf('加密完成 (%.4f 秒)。\n', encryptTime);
if ~isequal(sizeCheck, originalSize)
    warning('加密函数返回的尺寸不匹配！');
end

% --- 4. 获取用于分析/显示的 uint8 加密图像 ---
% 解码 DNA (不进行逆扩散) 以获得 uint8 表示
fprintf('生成用于分析的 uint8 加密图像...\n');
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
    fprintf('生成 uint8 加密图像时出错: %s\n', ME.message);
    return;
end
fprintf('Uint8 加密图像已生成。\n');

% --- 5. 显示图像 (原始, 加密) ---
figure('Name', '加密结果', 'NumberTitle', 'off'); % <--- 修改这里
subplot(1, 2, 1); imshow(originalImage); title('原始图像'); % <--- 修改这里
subplot(1, 2, 2); imshow(encryptedImageUint8); title('加密图像 (uint8)'); % <--- 修改这里

% --- 6. 直方图分析 ---
figure('Name', '直方图分析', 'NumberTitle', 'off'); % <-- (可选) 修改 Figure Name
if numChannels == 1
    subplot(1, 2, 1); imhist(originalImage); title('原始直方图'); ylim('auto'); % <-- (可选) 修改 title
    subplot(1, 2, 2); imhist(encryptedImageUint8); title('加密直方图'); ylim('auto'); % <-- (可选) 修改 title
else % Color: Show histograms for each channel
    colorLabels = {'红色', '绿色', '蓝色'}; % <-- (可选) 修改标签
    for k = 1:3
        subplot(2, 3, k);
        imhist(originalImage(:,:,k)); title(['原始', colorLabels{k}, '通道']); ylim('auto'); % <-- (可选) 修改 title
        subplot(2, 3, k+3);
        imhist(encryptedImageUint8(:,:,k)); title(['加密', colorLabels{k}, '通道']); ylim('auto'); % <-- (可选) 修改 title
    end
    % Adjust layout slightly
    annotation('textbox', [0 0.9 1 0.1], 'String', '直方图分析 (RGB 通道)', ... % <-- (可选) 修改 Annotation
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% --- 7. 信息熵分析 ---
fprintf('\n--- 信息熵分析 ---\n');
totalEntropyOrig = 0;
totalEntropyEnc = 0;
if numChannels == 1
    entropyOrig = calculateEntropy(originalImage);
    entropyEnc = calculateEntropy(encryptedImageUint8);
    fprintf('原始灰度熵: %.4f\n', entropyOrig);
    fprintf('加密灰度熵: %.4f (理想值 ~8)\n', entropyEnc);
    totalEntropyOrig = entropyOrig;
    totalEntropyEnc = entropyEnc;
else
    fprintf('通道熵:\n');
    fprintf('%10s %10s %10s\n', '通道', '原始', '加密');
    colorLabels = {'红', '绿', '蓝'}; % 使用上面定义的中文标签
    for k = 1:3
        entropyOrig = calculateEntropy(originalImage(:,:,k));
        entropyEnc = calculateEntropy(encryptedImageUint8(:,:,k));
        fprintf('%10s %10.4f %10.4f\n', colorLabels{k}, entropyOrig, entropyEnc);
        totalEntropyOrig = totalEntropyOrig + entropyOrig;
        totalEntropyEnc = totalEntropyEnc + entropyEnc;
    end
    fprintf('平均原始熵: %.4f\n', totalEntropyOrig / 3);
    fprintf('平均加密熵: %.4f (理想值 ~8)\n', totalEntropyEnc / 3);
end

% --- 8. 相关性分析 ---
fprintf('\n--- 相邻像素相关性分析 ---\n');
% 对彩色图像的每个通道进行分析
analyzeCorrelation(originalImage, '原始图像'); % <-- (可选) 修改标签
analyzeCorrelation(encryptedImageUint8, '加密图像'); % 分析 uint8 加密版本, <-- (可选) 修改标签

% --- 9. NPCR & UACI 分析 ---
fprintf('\n--- NPCR & UACI 分析 ---\n');
try
    % 创建一个轻微修改的原始图像 (翻转一个像素中的一个比特)
    originalImage_modified = originalImage;
    rand_row = randi([1, rows]);
    rand_col = randi([1, cols]);
    rand_channel = randi([1, numChannels]); % 选择一个通道进行修改
    originalImage_modified(rand_row, rand_col, rand_channel) = bitxor(originalImage(rand_row, rand_col, rand_channel), 1);

    % 使用相同的密钥加密修改后的图像
    fprintf('为 NPCR/UACI 加密修改后的图像...\n');
    [encryptedDnaData_modified, ~] = encryptImageDNA(originalImage_modified, keyStreamBinaryRule, keyStreamByteDiffusion);

    % 生成修改后加密图像的 uint8 版本
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
    fprintf('修改后的加密图像已生成。\n');

    % 计算 NPCR 和 UACI (函数内部处理颜色)
    [npcr_avg, uaci_avg] = calculateNPCR_UACI(encryptedImageUint8, encryptedImageUint8_modified);
    fprintf('平均 NPCR: %.4f%% (理想值 > 99.6%%)\n', npcr_avg);
    fprintf('平均 UACI: %.4f%% (对于 uint8 理想值 ~ 33.4%%)\n', uaci_avg);

catch ME
    fprintf('NPCR/UACI 计算过程中出错: %s\n', ME.message);
    fprintf('文件: %s, 行: %d\n', ME.stack(1).file, ME.stack(1).line);
end

% --- 10. 解密 ---
fprintf('\n解密图像...\n');
tic;
try
    % 传递加密的 DNA 数据 (cell 或 char) 和原始尺寸
    decryptedImage = decryptImageDNA(encryptedDnaData, originalSize, keyStreamBinaryRule, keyStreamByteDiffusion);
catch ME
    fprintf('解密过程中出错: %s\n', ME.message);
    fprintf('文件: %s, 行: %d\n', ME.stack(1).file, ME.stack(1).line);
    decryptedImage = []; % 标记解密失败
end
decryptTime = toc;

if ~isempty(decryptedImage)
    fprintf('解密完成 (%.4f 秒)。\n', decryptTime);

    % --- 11. 显示解密结果与验证 ---
    figure('Name', '解密结果与验证', 'NumberTitle', 'off'); % <--- 修改这里
    subplot(1, 2, 1); imshow(originalImage); title('原始图像'); % <--- 修改这里
    subplot(1, 2, 2); imshow(decryptedImage); title('解密图像'); % <--- 修改这里

    % 计算 MSE 和 PSNR
    fprintf('\n--- 解密质量评估 ---\n');
    try
        if numChannels == 1
            mse = immse(decryptedImage, originalImage);
            psnr_val = psnr(decryptedImage, originalImage);
        else % Color: 计算跨通道的平均 MSE/PSNR
            mse_r = immse(decryptedImage(:,:,1), originalImage(:,:,1));
            mse_g = immse(decryptedImage(:,:,2), originalImage(:,:,2));
            mse_b = immse(decryptedImage(:,:,3), originalImage(:,:,3));
            mse = mean([mse_r, mse_g, mse_b]);

            % PSNR 需要非零 MSE。处理完美重建的情况。
            if mse < 1e-10
                 psnr_val = Inf;
            else
                psnr_r = psnr(decryptedImage(:,:,1), originalImage(:,:,1));
                psnr_g = psnr(decryptedImage(:,:,2), originalImage(:,:,2));
                psnr_b = psnr(decryptedImage(:,:,3), originalImage(:,:,3));
                % 对 PSNR dB 值求平均是常见做法
                psnr_val = mean([psnr_r, psnr_g, psnr_b]);
            end
        end

        fprintf('均方误差 (MSE): %.4g\n', mse);
        if isinf(psnr_val) || psnr_val > 100 % 检查 Inf 或非常高的 PSNR
             fprintf('峰值信噪比 (PSNR): Inf dB\n');
             fprintf('成功: 解密图像与原始图像完全匹配！\n'); % <-- (可选) 修改文本
        else
            fprintf('峰值信噪比 (PSNR): %.4f dB\n', psnr_val);
            if mse > 1e-6 % 允许微小的浮点差异
                 fprintf('警告: 解密图像与原始图像略有不同！\n'); % <-- (可选) 修改文本
            else
                 fprintf('成功: 解密图像在视觉上与原始图像相同。\n'); % <-- (可选) 修改文本
            end
        end
    catch ME_psnr
        fprintf('无法计算 PSNR/MSE: %s\n', ME_psnr.message);
    end
else
    fprintf('解密失败。无法评估质量。\n');
end

% --- 12. 密钥空间分析 ---
fprintf('\n--- 密钥空间分析 ---\n');
analyzeKeySpace(initialConditions, [rho, sigma, beta]);

% --- 13. 密钥敏感性分析 ---
fprintf('\n--- 密钥敏感性分析 ---\n');
try
    % 定义分析所需的核心函数句柄
    encrypt_func_handle = @encryptImageDNA;
    decrypt_func_handle = @decryptImageDNA; % 这里不是严格需要，但良好实践
    keygen_func_handle = @generateLorenzKeyStream;
    dna_decode_func_handle = @dna_decode; % 需要获取 uint8 用于比较
    npcr_uaci_func_handle = @calculateNPCR_UACI;

    [npcr_sens, uaci_sens] = analyzeKeySensitivity(...
        originalImage, ...
        encrypt_func_handle, ...
        keygen_func_handle, ...
        dna_decode_func_handle, ...
        npcr_uaci_func_handle, ...
        initialConditions, rho, sigma, beta, ...
        key_sensitivity_delta);

    fprintf('对密钥变化的敏感性 (delta=%.0e):\n', key_sensitivity_delta);
    fprintf('  密文之间的 NPCR: %.4f%%\n', npcr_sens);
    fprintf('  密文之间的 UACI: %.4f%%\n', uaci_sens);
    if npcr_sens > 99 && uaci_sens > 30
        fprintf('  结果: 检测到高密钥敏感性 (良好！)。\n'); % <-- (可选) 修改文本
    else
        fprintf('  结果: 检测到低密钥敏感性 (潜在弱点！)。\n'); % <-- (可选) 修改文本
    end
catch ME_sens
    fprintf('密钥敏感性分析过程中出错: %s\n', ME_sens.message);
    fprintf('文件: %s, 行: %d\n', ME_sens.stack(1).file, ME_sens.stack(1).line);
end


% --- 14. 噪声攻击分析 ---
fprintf('\n--- 噪声攻击分析 ---\n');
if ~isempty(decryptedImage) % 仅当初始加密/解密成功时运行
    try
        % 根据类型选择噪声参数
        noise_type_cn = ''; % 用于标题的中文名称
        if strcmpi(noise_type, 'gaussian')
            noise_param = noise_param_gaussian;
            noise_type_cn = '高斯';
            fprintf('应用高斯噪声 (方差 = %.4f)...\n', noise_param);
        elseif strcmpi(noise_type, 'salt & pepper')
            noise_param = noise_param_sp;
             noise_type_cn = '椒盐';
            fprintf('应用椒盐噪声 (密度 = %.2f)...\n', noise_param);
        else
            error('不支持的噪声类型: %s', noise_type);
        end

        % 将必要组件传递给分析函数
        [decryptedNoisyImage, psnr_noise, mse_noise] = analyzeNoiseAttack(...
            originalImage, ...
            encryptedImageUint8, ... % 攻击 uint8 版本
            keyStreamBinaryRule, ...
            keyStreamByteDiffusion, ...
            noise_type, noise_param, ...
            @dna_encode, ... % 传递函数句柄
            @decryptImageDNA); % 传递函数句柄

        fprintf('噪声攻击分析完成。\n');
        fprintf('  MSE (原始 vs. 解密的含噪图像): %.4f\n', mse_noise);
        fprintf('  PSNR (原始 vs. 解密的含噪图像): %.4f dB\n', psnr_noise);

        figure('Name', '噪声攻击结果', 'NumberTitle', 'off'); % <-- (可选) 修改 Figure Name
        subplot(1, 2, 1); imshow(originalImage); title('原始图像'); % <-- 使用中文标题
        subplot(1, 2, 2); imshow(decryptedNoisyImage); title(sprintf('%s噪声攻击后解密', noise_type_cn)); % <--- 修改这里

    catch ME_noise
        fprintf('噪声攻击分析过程中出错: %s\n', ME_noise.message);
         fprintf('文件: %s, 行: %d\n', ME_noise.stack(1).file, ME_noise.stack(1).line);
    end
else
    fprintf('由于初始解密失败，跳过噪声攻击分析。\n');
end


% --- 15. 裁剪攻击分析 ---
fprintf('\n--- 裁剪攻击分析 ---\n');
if ~isempty(decryptedImage) % 仅当初始加密/解密成功时运行
    try
        fprintf('应用裁剪攻击 (从中心裁剪 %.1f%%)...\n', crop_percentage*100);

        % 定义裁剪矩形 (居中)
        crop_width = floor(cols * crop_percentage);
        crop_height = floor(rows * crop_percentage);
        start_x = floor((cols - crop_width) / 2) + 1;
        start_y = floor((rows - crop_height) / 2) + 1;
        % 将坐标限制在图像边界内
        start_x = max(1, start_x);
        start_y = max(1, start_y);
        end_x = min(cols, start_x + crop_width - 1);
        end_y = min(rows, start_y + crop_height - 1);

        cropRectangle = [start_x, start_y, end_x - start_x + 1, end_y - start_y + 1]; % 某些函数期望的 [x, y, width, height] 格式

        % 将必要组件传递给分析函数
        [decryptedCroppedImage, psnr_crop, mse_crop] = analyzeCroppingAttack(...
            originalImage, ...
            encryptedImageUint8, ... % 攻击 uint8 版本
            keyStreamBinaryRule, ...
            keyStreamByteDiffusion, ...
            cropRectangle, ... % 传递计算出的矩形
            @dna_encode, ... % 传递函数句柄
            @decryptImageDNA); % 传递函数句柄

         fprintf('裁剪攻击分析完成。\n');
         fprintf('  MSE (原始 vs. 解密的裁剪图像): %.4f\n', mse_crop);
         fprintf('  PSNR (原始 vs. 解密的裁剪图像): %.4f dB\n', psnr_crop);

        figure('Name', '裁剪攻击结果', 'NumberTitle', 'off'); % <-- (可选) 修改 Figure Name
        subplot(1, 2, 1); imshow(originalImage); title('原始图像'); % <-- 使用中文标题
        subplot(1, 2, 2); imshow(decryptedCroppedImage); title(sprintf('%.0f%% 裁剪攻击后解密', crop_percentage*100)); % <--- 修改这里

    catch ME_crop
        fprintf('裁剪攻击分析过程中出错: %s\n', ME_crop.message);
        fprintf('文件: %s, 行: %d\n', ME_crop.stack(1).file, ME_crop.stack(1).line);
    end
else
     fprintf('由于初始解密失败，跳过裁剪攻击分析。\n');
end

fprintf('\n--- 所有处理和分析完成 ---\n');
