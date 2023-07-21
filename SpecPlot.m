function SpecPlot(signal,fs,bandtype)
% signal为输入信号
% fs为采样频率(高于两倍信号频率)
% bandtype为频谱绘制选择信号，'single'为单边频谱，'double'为双边频谱
if nargin < 2
    fprintf('error,使用open SpecPlot命令语句进行查看\n');
    return
end
if nargin == 2 % 默认画双边频谱
    bandtype = 'double';
end
len = length(signal); % 信号长度
mag = abs(fftshift(fft(signal))); % 频谱幅度|F(e^jw)|
mag = mag/max(mag); % 归一化
if (strcmp(bandtype,'single'))
    f = [0:len/2]*fs/len; % 所画频谱横坐标范围为[0:fs/2]
    figure
    semilogy(f,mag(round(len/2):len));title('频谱图');xlabel('f/Hz');ylabel('|X(f)|');
else
    f = [-len/2:len/2-1]*fs/len; % 所画频谱横坐标范围为[-fs/2:fs/2]
    figure
    semilogy(f,mag);title('频谱图');xlabel('f/Hz');ylabel('|X(f)|');
end
end
