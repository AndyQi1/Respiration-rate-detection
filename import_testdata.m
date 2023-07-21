%% 导入数据
Y = ['1'];
X = {'001','002','003','004','005','006','007','008','009'};
RMSE = zeros(1,length(X));
RMSE_baseline = zeros(1,length(X));
for SampleNum = 1%:length(X)
    name = strcat('TestCfgData',X{SampleNum},'.txt'); % 实例参数
    Cfgdata = importdata(strcat('TestData/',name));
    %     name = strcat('TestInputData',X{SampleNum},'.txt'); % CSI数据
    load(strcat('TestInputData',X{SampleNum},'.mat'))
    %     Inputdata = importdata(strcat('TestData/',name));
    name = strcat('TestBreathWave',X{SampleNum},'.txt'); % 呼吸波形
    Breathwave = importdata(strcat('TestData/',name));
    name = strcat('TestGroundTruthData',X{SampleNum},'.txt'); % 呼吸波形
    GroundTruth = importdata(strcat('TestData/',name));
    
    % Inputdata = reshape(Inputdata', 1, []); % 转换为一个行向量
    dataNum = Cfgdata(1); % 数据数
    dataPerson = Cfgdata(2:dataNum + 1); % 每条数据的人数
    N_Tx = Cfgdata(dataNum + 2); % 发射天线数
    N_Rx = Cfgdata(dataNum + 3); % 接收天线数
    N_Sc = Cfgdata(dataNum + 4); % 子载波数
    N_T = Cfgdata(dataNum + 5:2*dataNum + 4); % 测量次数
    T_Dur = Cfgdata(2*dataNum + 5:3*dataNum + 4); % 采集持续时间
    f_Start = Cfgdata(end - 1); % 起始频率
    f_End = Cfgdata(end); % 终止频率
    fs = (N_T - 1)./T_Dur; % 采样频率
    f_Center = (f_Start + f_End)/2; % 中心频率
    delta_f = (f_End - f_Start)/N_Sc;
    % 数据每列先遍历接收天线，再遍历子载波， 每行代表时序上的信息
    Idx = cumsum([0;N_T]);
    BPM = zeros(dataNum,3);
    
    for ii = 1:length(N_T)
        
        data = Inputdata(Idx(ii) + 1:Idx(ii + 1),:);
        real = data(:,1:2:end);
        imag = data(:,2:2:end);
        data = real + 1j*imag;
       
        
        phase = zeros(N_T(ii),N_Sc*(N_Rx - 1));
        for ss = 1:N_Sc
            for nn = 2:N_Rx
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = 180*angle(data(:,(ss-1)*N_Rx + nn))/pi - 180*angle(data(:,(ss-1)*N_Rx + 1))/pi; % 转为角度制
            end
        end
        Type = 'db4';
        level = 4;
        close all;
        col_num = 2;
        SpecPlot(phase(:,col_num),fs(ii),'double');
        PhaseTrend = hampel(phase(:,1),2000,0.01);
        phase(:,col_num) = phase(:,col_num) - PhaseTrend;
        SpecPlot(phase(:,col_num),fs(ii),'double');
        phase(:,col_num) = hampel(phase(:,col_num));
        SpecPlot(phase(:,col_num),fs(ii),'double');
        for kk = 1:N_Sc*(N_Rx - 1)
            phase(:,kk) = hampel(phase(:,kk));
            [c, l] = wavedec(phase(:,kk), level, Type);
            c(1:l(3)) = 0;
            phase(:,kk) = waverec(c,l,'db4');
        end
        % 子载波筛选
        ChosenNum = 3;
        MAD = mad(phase,0,1);
        [~,MADInd] = sort(MAD,'descend');
        ChosenPhase = sum(phase(:,MADInd(1:end)));
        
        phaseFre = abs(fftshift(fft(ChosenPhase)));
        len = length(phaseFre); % 相位数据长度
        
        [~,Index] = sort(phaseFre(round(len/2):end),'descend');
        LowerIndex = round(5*len/fs(ii)); % 呼吸率在5-50Hz之间，选取该范围的频率最大值
        UpperIndex = round(50*len/fs(ii));
        Index = Index(Index >= LowerIndex);
        Index = Index(Index <= UpperIndex);
        BPM(ii,1:dataPerson(ii)) = Index(1:dataPerson(ii))'/len*fs(ii);
        ii
        %     close all;
        %     figure(1)
        %     plot(1:length(data),phaseFre);
        %     SpecPlot(phase, fs(ii), 'single')
    end
    
    
    BPM = reshape(BPM',1,[]);
    BPM = BPM(BPM > 0);
    GroundTruth(isnan(GroundTruth)) = 0;
    GroundTruth = reshape(GroundTruth',1,[]);
    GroundTruth = GroundTruth(GroundTruth > 0);
    %
    RMSE(SampleNum) = sqrt(1/dataNum*sum((GroundTruth - BPM).^2));
    %% 基准
    BPM_baseline = zeros(dataNum,3);
    for ii = 1:length(N_T)
        
        data = Inputdata(Idx(ii) + 1:Idx(ii + 1),:);
        real = data(:,1:2:end);
        imag = data(:,2:2:end);
        data = real + 1j*imag;
        phase_baseline = zeros(N_T(ii),1);
        for ss = 1:N_Sc
            for nn = 2:N_Rx
                phase_baseline = 180*angle(data(:,(ss-1)*N_Rx + nn))/pi - 180*angle(data(:,(ss-1)*N_Rx + 1))/pi + phase_baseline; % 转为角度制
            end
        end
        phase_baseline = 1/(N_Sc*N_Rx)*phase_baseline; % 归一化
        
        phase_baseline = hampel(phase_baseline);
        %         len = length(data); % 数据长度
        % %         小波变换
        
        [c, l] = wavedec(phase_baseline, 4, 'db4');
        %         phase = appcoef(c ,l, 'db4', 4);
        coef = detcoef(c ,l, 'db4', 4);
        phase_baseline = coef{3};
        %         phase = hampel(phase);
        
        phaseFre = abs(fftshift(fft(phase_baseline)));
        len = length(phaseFre); % 相位数据长度
        
        [~,Index] = sort(phaseFre(round(len/2):end),'descend');
        LowerIndex = round(5*len/fs(ii)); % 呼吸率在5-50Hz之间，选取该范围的频率最大值
        UpperIndex = round(50*len/fs(ii));
        Index = Index(Index >= LowerIndex);
        Index = Index(Index <= UpperIndex);
        BPM_baseline(ii,1:dataPerson(ii)) = Index(1:dataPerson(ii))'/len*fs(ii);
        
        %     close all;
        %     figure(1)
        %     plot(1:length(data),phaseFre);
        %     SpecPlot(phase, fs(ii), 'single')
    end
    
    
    BPM_baseline = reshape(BPM_baseline',1,[]);
    BPM_baseline = BPM_baseline(BPM_baseline > 0);
    
    %
    RMSE_baseline(SampleNum) = sqrt(1/dataNum*sum((GroundTruth - BPM_baseline).^2));
    SampleNum
end
Type
RMSE_baseline
RMSE
