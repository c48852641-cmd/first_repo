%% ============================================================
%  TDMS信号分析仪 (全血优化·UI交互·ADC物理功率映射版)
%  适用：
%    - 各种宽带连续信号与短促突发信号 (LTE, 5G NR, Wi-Fi, 蓝牙等)
%    - 复杂多载波/强干扰环境的可视化排查
%    - 处理 GB 级别的海量 TDMS 原始量化采集数据
%  核心算法升级：
%    - 引入 Peak-Hold STFT 算法防 Burst 信号稀释
%    - 严格区分绝对功率 (dBm) 与 功率谱密度 (dBm/Hz)
% ============================================================

clc; clear; close all;

%% 1. 用户配置 (UI 交互版)

% --- 1.1 弹出文件选择框 ---
[file_name, folder_path] = uigetfile('*.tdms', '请选择要分析的 TDMS 数据文件');
if isequal(file_name, 0)
    disp('用户取消了文件选择，程序退出。');
    return; % 安全退出
end
file_path = fullfile(folder_path, file_name);

% --- 1.2 弹出参数输入对话框 ---
prompt = {
    '基带采样率 Fs (单位: MS/s):', ...
    '中心频率 Fc (单位: MHz):', ...
    '参考电平 Ref Level (单位: dBm):'
};
dlgtitle = 'TDMS 信号分析仪 - 参数设置';
dims = [1 55];

% 默认值
definput = {'50', '2460', '-40'};

% 呼出对话框
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('用户取消了参数输入，程序退出。');
    return;
end

% --- 1.3 解析用户输入 ---
Fs = str2double(answer{1}) * 1e6;     
Fc = str2double(answer{2}) * 1e6;     
ref_level_dbm = str2double(answer{3}); 

% --- 1.4 高级系统配置 (ADC 硬件物理参数) ---
adc_bit_depth = 16;                        % NI-RFSA 驱动存储位宽为 16-bit
adc_full_scale = 2^(adc_bit_depth - 1);    % 对应满量程值 32768 (0 dBFS)

chunk_size = 5e6;          % 分块读取大小
MAX_FILE_SIZE_GB = 2;      % 最大安全文件大小
MAX_ESTIMATED_RAM_GB = 12.0;

nfft_psd = 8192;           % 整体频谱分辨率
nfft_spec = 512;           % 时频图分辨率 (兼顾时间与频率)
num_time_bins = 800;       % 时频图时间切片数

%% 2. 文件与内存检查
finfo = dir(file_path);
file_size_bytes = finfo.bytes;
file_size_mb = file_size_bytes / 1024^2;
file_size_gb = file_size_bytes / 1024^3;

fprintf('\n=========== 1. 文件检查 ===========\n');
fprintf('文件路径: %s\n', file_path);
fprintf('文件大小: %.3f MB (%.6f GB)\n', file_size_mb, file_size_gb);
fprintf('采样率: %g MS/s\n', Fs/1e6);
fprintf('中心频率: %g MHz\n', Fc/1e6);
fprintf('参考电平: %g dBm\n', ref_level_dbm);

if file_size_gb > MAX_FILE_SIZE_GB
    error('文件大小 %.3f GB 超过阈值 %.3f GB，程序退出以防死机。', file_size_gb, MAX_FILE_SIZE_GB);
end

estimated_ram_gb = file_size_gb * 12;
fprintf('\n预计 MATLAB 峰值内存占用: %.3f GB\n', estimated_ram_gb);

if exist('memory','file') == 2
    [usr, ~] = memory;
    avail_gb = usr.MemAvailableAllArrays / 1024^3;
    fprintf('当前 MATLAB 可用内存: %.3f GB\n', avail_gb);
    if estimated_ram_gb > 0.7 * avail_gb
        error('预计内存占用接近当前可用内存，程序退出。');
    end
end

% 初始化进度条
h_wait = waitbar(0, '正在初始化...', 'Name', 'TDMS 数据处理进度');

%% 3. 极速读取 IQ 数据与预处理
fprintf('\n=========== 2. 读取 TDMS 数据 ===========\n');
waitbar(0.05, h_wait, '正在读取 TDMS 数据 (分块加载)...');

ds = tdmsDatastore(file_path);
ds.ReadSize = chunk_size;

chunk_list = cell(1000, 1); 
chunk_count = 0;

while hasdata(ds)
    T = read(ds);
    if iscell(T)
        T = T{1};
    end
    raw_vec = double(T.Variables);
    if size(raw_vec,2) > 1
        raw_vec = raw_vec(:,1);
    end
    I = raw_vec(1:2:end);
    Q = raw_vec(2:2:end);
    len = min(length(I), length(Q));
    if len == 0, continue; end
    
    chunk_count = chunk_count + 1;
    if chunk_count > length(chunk_list)
        chunk_list = [chunk_list; cell(500, 1)]; %#ok<AGROW> 
    end
    chunk_list{chunk_count} = complex(I(1:len), Q(1:len)); 
    
    if mod(chunk_count, 5) == 0 && ishandle(h_wait)
        waitbar(0.05 + 0.25 * min(1, chunk_count/200), h_wait, sprintf('已读取 %d 个数据块...', chunk_count));
    end
end

if chunk_count == 0
    close(h_wait);
    error('IQ 数据读取失败，结果为空。');
end

waitbar(0.3, h_wait, '正在拼接数据与去除直流偏置...');
sig_all = vertcat(chunk_list{1:chunk_count});
clear chunk_list; 

% 只消除绝对 0Hz 的静态直流偏置
sig_clean = sig_all - mean(sig_all);
L = length(sig_clean);

fprintf('复采样点数: %d\n', L);
fprintf('采样时长: %.6f s\n', L/Fs);

%% 4. 整体频谱计算与 ADC 诊断 (修正为真 PSD 计算)
fprintf('\n=========== 3. 信号全景计算 ===========\n');
if ishandle(h_wait), waitbar(0.4, h_wait, '正在计算全局功率谱...'); end

peak_amplitude = max(abs(sig_clean));
fprintf('当前信号原始量化峰值: %.1f (硬件满量程边界: %d)\n', peak_amplitude, adc_full_scale);
if peak_amplitude > adc_full_scale * 0.99
    warning('警告: 检测到 ADC 可能发生削波 (Clipping)！请在收发器调高 Reference Level。');
end

% 满量程归一化
sig_psd_est = sig_clean(1:min(L, 5*Fs)) / adc_full_scale; 

% 严格使用 'psd' 计算功率谱密度 (单位 1/Hz)
[pxx_psd, f] = pwelch(sig_psd_est, hann(nfft_psd), floor(nfft_psd/2), nfft_psd, Fs, 'centered', 'psd');

pxx_psd = suppress_center_spur(pxx_psd, 6); 
% 相对满量程对数运算后加上 Reference Level，得到纯正的 dBm/Hz
pxx_db_absolute = 10*log10(pxx_psd + eps) + ref_level_dbm;

%% 5. 时频图计算 (修正为 Peak-Hold STFT)
fprintf('\n=========== 4. 时频图计算 ===========\n');
if ishandle(h_wait), waitbar(0.5, h_wait, '正在计算时频矩阵...'); end

step_size = max(1, floor(L / num_time_bins));
spec_matrix = zeros(nfft_spec, num_time_bins);

for k = 1:num_time_bins
    idx1 = (k-1)*step_size + 1;
    idx2 = min(k*step_size, L);
    
    % 切片数据除以满量程进行归一化
    segment = sig_clean(idx1:idx2) / adc_full_scale;

    if length(segment) < 64
        continue;
    end

    % 使用 spectrogram 计算该时间块内所有的短时频谱快照
    seg_win = min(nfft_spec, length(segment));
    seg_ov  = floor(seg_win / 2);
    [~, ~, ~, P_matrix] = spectrogram(segment, hann(seg_win), seg_ov, nfft_spec, Fs, 'centered', 'power');
    
    % 核心算法升级：峰值保持 (Max Hold)，提取每一频点在该时间块内的最大能量
    P_peak_hold = max(P_matrix, [], 2);
    
    P_peak_hold = suppress_center_spur(P_peak_hold, 5); 
    spec_matrix(:,k) = 10 * log10(P_peak_hold + eps);
    
    if mod(k, 40) == 0 && ishandle(h_wait)
        progress = 0.5 + 0.4 * (k / num_time_bins);
        waitbar(progress, h_wait, sprintf('时频图处理中... %d%%', round((k/num_time_bins)*100)));
    end
end

f_axis = linspace(-Fs/2, Fs/2, nfft_spec)/1e6 + Fc/1e6;
t_axis = linspace(0, L/Fs, num_time_bins);

%% 6. 可视化绘图
fprintf('\n=========== 5. 绘图输出 ===========\n');
if ishandle(h_wait), waitbar(0.95, h_wait, '正在渲染可视化图表...'); end

fig = figure('Color','w','Position',[80 50 1400 850]);
sgtitle(sprintf('宽带信号基础分析\n%s', file_name), 'Interpreter', 'none', 'FontWeight', 'bold');

% --- 1. 宽带时域包络 (映射到绝对功率 dBm) ---
subplot(2,2,1);
decim = max(1, floor(L / 200000)); 
sig_env_raw = abs(sig_clean(1:decim:end)); 
t_env = (0:length(sig_env_raw)-1) * decim / Fs;

sig_env_norm = sig_env_raw / adc_full_scale;
sig_env_dbm = ref_level_dbm + 20*log10(sig_env_norm + eps);

plot(t_env, sig_env_dbm, 'Color', [0 0.4470 0.7410]);
grid on; axis tight;
title(sprintf('1. 宽带时域功率包络 (Ref Level: %g dBm)', ref_level_dbm));
xlabel('时间 (s)');
ylabel('瞬时绝对功率 (dBm)');
ylim([max(min(sig_env_dbm), ref_level_dbm - 80), ref_level_dbm + 5]); 

% --- 2. 整体绝对功率谱密度 (dBm/Hz) ---
subplot(2,2,2);
plot((f + Fc)/1e6, pxx_db_absolute, 'r', 'LineWidth', 1);
grid on; axis tight;
title('2. 整体功率谱密度 (真实底噪评估)');
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dBm/Hz)'); % 单位已修正
xlim([Fc-Fs/2, Fc+Fs/2]/1e6);

% --- 3. 时频瀑布图 (峰值保持版) ---
subplot(2,1,2);
% 加上参考电平偏置，显示突发信号的峰值绝对功率
spec_matrix_absolute = spec_matrix + ref_level_dbm;
imagesc(t_axis, f_axis, spec_matrix_absolute);
axis xy; colormap(jet); colorbar;
title('3. 时频瀑布图 (Peak-Hold 峰值保持，绝对功率 dBm/Bin)');
xlabel('时间 (s)');
ylabel('频率 (MHz)');

vals = spec_matrix_absolute(:);
caxis([prctile(vals, 50), prctile(vals, 99.8)]); 

if ishandle(h_wait)
    waitbar(1.0, h_wait, '处理完成！');
    pause(0.5);
    close(h_wait);
end
fprintf('\n基础分析渲染完成。\n');

% 渲染完成后，打印详细元数据
print_tdms_metadata(file_path);

%% ======================== 本地函数 ========================

function pxx_out = suppress_center_spur(pxx_in, spur_bins)
    % 消除中心频点的本振泄漏伪影
    pxx_out = pxx_in;
    N = length(pxx_in);
    c = floor(N/2) + 1;

    idx1 = max(1, c - spur_bins);
    idx2 = min(N, c + spur_bins);

    left_idx  = max(1, idx1 - 5);
    right_idx = min(N, idx2 + 5);

    left_vals  = pxx_in(left_idx:idx1-1);
    right_vals = pxx_in(idx2+1:right_idx);

    neigh = [left_vals(:); right_vals(:)];
    if isempty(neigh)
        return;
    end

    fill_val = median(neigh);
    pxx_out(idx1:idx2) = fill_val;
end

function print_tdms_metadata(file_path)
    % 打印 TDMS 元数据（保持原样，极其健壮）
    fprintf('\n================ TDMS 元数据完整打印 ================\n');

    if ~isfile(file_path)
        fprintf('文件不存在: %s\n', file_path);
        return;
    end

    finfo = dir(file_path);
    fprintf('[文件基本信息]\n');
    fprintf('文件路径: %s\n', file_path);
    fprintf('文件名称: %s\n', finfo.name);
    fprintf('文件大小: %.3f MB (%.6f GB)\n', finfo.bytes/1024^2, finfo.bytes/1024^3);
    fprintf('修改时间: %s\n', finfo.date);

    try
        info = tdmsinfo(file_path);
    catch ME
        fprintf('tdmsinfo 读取失败: %s\n', ME.message);
        return;
    end

    fprintf('\n[tdmsinfo 顶层字段]\n');
    try
        top_props = properties(info);
    catch
        if isstruct(info)
            top_props = fieldnames(info);
        else
            fprintf('无法解析 tdmsinfo 顶层字段。\n');
            disp(info);
            return;
        end
    end

    for i = 1:length(top_props)
        name = top_props{i};
        try
            value = info.(name);
            print_any_value(name, value, 2);
        catch
            fprintf('  %s : [读取失败]\n', name);
        end
    end

    fprintf('\n[ChannelList 详细信息]\n');
    try
        if isprop(info, 'ChannelList') || (isstruct(info) && isfield(info, 'ChannelList'))
            chList = info.ChannelList;
            if istable(chList)
                fprintf('ChannelList 大小: %d x %d\n', size(chList,1), size(chList,2));
                disp(chList);

                fprintf('\n[逐通道详细信息]\n');
                varNames = chList.Properties.VariableNames;
                for r = 1:height(chList)
                    fprintf('--------------------------------------------------\n');
                    fprintf('通道 #%d\n', r);
                    for c = 1:length(varNames)
                        vn = varNames{c};
                        val = chList{r,c};
                        fprintf('  %s : ', vn);
                        print_scalar_or_brief(val);
                    end
                end
            else
                fprintf('ChannelList 不是 table，内容如下：\n');
                disp(chList);
            end
        else
            fprintf('未发现 ChannelList 字段。\n');
        end
    catch ME
        fprintf('打印 ChannelList 失败: %s\n', ME.message);
    end

    fprintf('\n[PropertyList 详细信息]\n');
    try
        if isprop(info, 'PropertyList') || (isstruct(info) && isfield(info, 'PropertyList'))
            propList = info.PropertyList;
            if istable(propList)
                fprintf('PropertyList 大小: %d x %d\n', size(propList,1), size(propList,2));
                disp(propList);

                fprintf('\n[逐属性打印]\n');
                varNames = propList.Properties.VariableNames;
                for r = 1:height(propList)
                    fprintf('--------------------------------------------------\n');
                    fprintf('属性 #%d\n', r);
                    for c = 1:length(varNames)
                        vn = varNames{c};
                        val = propList{r,c};
                        fprintf('  %s : ', vn);
                        print_scalar_or_brief(val);
                    end
                end
            else
                fprintf('PropertyList 不是 table，内容如下：\n');
                disp(propList);
            end
        else
            fprintf('未发现 PropertyList 字段。\n');
        end
    catch ME
        fprintf('打印 PropertyList 失败: %s\n', ME.message);
    end

    fprintf('\n[按通道组(Group)整理输出]\n');
    try
        if exist('chList','var') && istable(chList)
            vns = lower(chList.Properties.VariableNames);

            idxGroup = find(contains(vns,'group'), 1);
            idxChannel = find(contains(vns,'channel'), 1);
            idxNumSamples = find(contains(vns,'numsamples') | contains(vns,'sample'), 1);
            idxDataType = find(contains(vns,'datatype') | contains(vns,'type'), 1);

            if isempty(idxGroup) || isempty(idxChannel)
                fprintf('ChannelList 中未能自动识别 Group/Channel 列。\n');
            else
                groups = chList{:, idxGroup};
                channels = chList{:, idxChannel};

                if iscell(groups)
                    ugroups = unique(groups);
                else
                    ugroups = unique(groups, 'stable');
                end

                for g = 1:numel(ugroups)
                    fprintf('==================================================\n');
                    fprintf('通道组: ');
                    print_scalar_or_brief(ugroups(g));

                    mask = false(height(chList),1);
                    for rr = 1:height(chList)
                        if iscell(groups)
                            mask(rr) = isequaln(groups{rr}, ugroups{g});
                        else
                            mask(rr) = isequaln(groups(rr), ugroups(g));
                        end
                    end

                    rows = find(mask);
                    for rr = rows(:).'
                        fprintf('  通道: ');
                        if iscell(channels)
                            print_scalar_or_brief(channels{rr});
                        else
                            print_scalar_or_brief(channels(rr));
                        end

                        if ~isempty(idxNumSamples)
                            fprintf('    NumSamples: ');
                            print_scalar_or_brief(chList{rr, idxNumSamples});
                        end
                        if ~isempty(idxDataType)
                            fprintf('    DataType: ');
                            print_scalar_or_brief(chList{rr, idxDataType});
                        end
                    end
                end
            end
        else
            fprintf('没有可用于分组整理的 ChannelList。\n');
        end
    catch ME
        fprintf('按通道组整理输出失败: %s\n', ME.message);
    end

    fprintf('\n================ 元数据打印结束 ================\n');
end

function print_any_value(name, value, indent)
    sp = repmat(' ',1,indent);

    if istable(value)
        fprintf('%s%s : [table %d x %d]\n', sp, name, size(value,1), size(value,2));
    elseif isstruct(value)
        fprintf('%s%s : [struct]\n', sp, name);
        fn = fieldnames(value);
        for k = 1:length(fn)
            try
                print_any_value(fn{k}, value.(fn{k}), indent+2);
            catch
                fprintf('%s  %s : [读取失败]\n', sp, fn{k});
            end
        end
    elseif iscell(value)
        fprintf('%s%s : [cell %d x %d]\n', sp, name, size(value,1), size(value,2));
    elseif isnumeric(value) || islogical(value)
        if isscalar(value)
            fprintf('%s%s : %g\n', sp, name, value);
        else
            fprintf('%s%s : [%s %s]\n', sp, name, class(value), mat2str(size(value)));
        end
    elseif ischar(value) || isstring(value)
        fprintf('%s%s : %s\n', sp, name, string(value));
    else
        fprintf('%s%s : [%s]\n', sp, name, class(value));
    end
end

function print_scalar_or_brief(val)
    if iscell(val)
        if numel(val) == 1
            print_scalar_or_brief(val{1});
        else
            fprintf('[cell %d elements]\n', numel(val));
        end
        return;
    end

    if isstring(val) || ischar(val)
        fprintf('%s\n', string(val));
    elseif isnumeric(val) || islogical(val)
        if isscalar(val)
            fprintf('%g\n', val);
        else
            fprintf('[%s %s]\n', class(val), mat2str(size(val)));
        end
    else
        try
            disp(val);
        catch
            fprintf('[%s]\n', class(val));
        end
    end
end