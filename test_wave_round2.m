%% 导入数据
Y = ['1'];
X = {'001','002','003','004','005','006','007','008','009'};
RMSE = zeros(1,length(X));
RMSE_baseline = zeros(1,length(X));
MinBPM = 5;
MaxBPM = 45;
DataNum = 0;
RMSE_sum = 0;
close all;
plotfig = 0;
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
    Type = 'db4';
    level = 4;
    for ii = 1%:length(N_T)
        
        data = Inputdata(Idx(ii) + 1:Idx(ii + 1),:);
        real = data(:,1:2:end);
        imag = data(:,2:2:end);
        data = real + 1j*imag;
        pNum = N_Sc*N_Rx*(N_Rx - 1)/2;
        phase = zeros(N_T(ii),pNum);
        
        len = N_T(ii);
        WindowLen = round(len/(T_Dur(ii)*(MaxBPM + MinBPM)/120)); % 计算滑动窗口
        LoadPhase = cell(1,pNum);
        pr = zeros(1,pNum);
        for ss = 1:N_Sc
            for ff = 1:N_Rx - 1
            for nn = ff + 1:N_Rx
                fff = ff - 1;
                phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff) =...
                    180*angle(conj(data(:,(ss-1)*N_Rx + ff)).*data(:,(ss-1)*N_Rx + nn))/pi; % 转为角度制
                phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff) = detrend(phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff), 1); % 去除趋势
                phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff) = hampel(phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff), 10);
        
                Movestd = movstd(phase(:,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff),WindowLen);
                Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
                StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
                EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
                prefix = 50;

%                 LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = phase([StartPoint + Precyclix:EndPoint - Precyclix]...
%                     ,(ss-1)*(N_Rx-1) + nn - 1);
                
                LoadPhase{(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff} = normalize(phase([StartPoint + prefix:EndPoint - prefix]...
                    ,(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff));
                if plotfig == 1
                    figure(1);plot(LoadPhase{(ss-1)*(N_Rx*(N_Rx - 1)/2) + (2*N_Rx - fff - 1)*fff/2 + nn - ff});hold on;
                end
            end
            end
            
        end
        
        % 子载波筛选
%         [~,ChosenIndex] = max(pr);
        ChosenStart = 1;
        ChosenNum = 20;
        MAD = zeros(1, pNum);
        for mm = 1:pNum
            MAD(mm) = mad(LoadPhase{mm},0);
        end
        fpass = [MinBPM, MaxBPM]/60; % 通带频率范围
        % % % % 滤波信号
        
        [~,MADInd] = sort(MAD,'descend');
        pkssum = zeros(1,ChosenNum);
        for cc = ChosenStart:ChosenStart + ChosenNum - 1
            LoadPhase{MADInd(cc)} = bandpass(LoadPhase{MADInd(cc)}, fpass, fs(ii));
            Phasecorr = xcorr(LoadPhase{MADInd(cc)});
            [~,Slefcorr] = max(Phasecorr);
            Phasecorr = Phasecorr([1:Slefcorr - 1,Slefcorr + 1:length(Slefcorr)]);
            Phasecorr = Phasecorr/max(Phasecorr);
            [pks, ~, ~, ~] = findpeaks(Phasecorr);
            pks = sort(pks,'descend');
            pkssum(cc) = sum(pks(1:min(length(pks),10)));
        end
        [~,pkssumIndex] = sort(pkssum,'descend');
        
        % 呼吸率估算
        ChosenIndex = MADInd(pkssumIndex(1:ChosenNum)); 
        SelPhase = [];
        SelNum = 5; % 选取三条
        WaveLen = 0.45;
        
        if dataPerson(ii) == 1
            for ss = 1:SelNum
                SelPhase = LoadPhase{ChosenIndex(ss)};
                [timepeak, Index, w_phase, p_phase] = findpeaks(SelPhase); % 时域上峰值
                w_phase = w_phase(timepeak >= mean(SelPhase) + 0.05*sqrt(var(SelPhase)));
                p_phase = p_phase(timepeak >= mean(SelPhase) + 0.05*sqrt(var(SelPhase)));
                Index = Index(timepeak >= mean(SelPhase) + 0.05*sqrt(var(SelPhase)));
                Index = Index(w_phase >= WaveLen*fs(ii)*60/MaxBPM);
                timepeakNum = sum(- Index(1:end-1) + Index(2:end) >= 60/MaxBPM*fs(ii))  + 1;
                BPM_per_Selphase(ss) = 60/(length(SelPhase)/fs(ii)/timepeakNum);
                len = length(SelPhase);
                SelPhase = reshape(SelPhase,[],1);
                x_axis = [-len/2:len/2-1]*fs(ii)*60/len;
                if plotfig == 1
                    figure;subplot 211;
                    plot(SelPhase); hold on; plot(Index,SelPhase(Index),'r*');hold on;
                    subplot 212; plot(x_axis,abs(fftshift(fft(SelPhase))));
                end
            end
            BPM_in_Time = mean(BPM_per_Selphase);
            BPM(ii,1:dataPerson(ii)) = BPM_in_Time;
        else
            for ss = 1:SelNum
                SelPhase1 = LoadPhase{ChosenIndex(ss)};
                [timepeak, Index, w_phase, p_phase] = findpeaks(SelPhase1); % 时域上峰值
                w_phase = w_phase(timepeak >= mean(SelPhase1) + 0.05*sqrt(var(SelPhase1)));
                p_phase = p_phase(timepeak >= mean(SelPhase1) + 0.05*sqrt(var(SelPhase1)));
                Index = Index(timepeak >= mean(SelPhase1) + 0.05*sqrt(var(SelPhase1)));
                Index = Index(w_phase >= WaveLen*fs(ii)*60/MaxBPM);
                timepeakNum = sum(- Index(1:end-1) + Index(2:end) >= 60/MaxBPM*fs(ii))  + 1;
                BPM_per_Selphase(ss) = 60/(length(SelPhase1)/fs(ii)/timepeakNum);
                SelPhase1 = reshape(SelPhase1,[],1);
                
                len = length(SelPhase1);
                x_axis = [-len/2:len/2-1]*fs(ii)*60/len;
                if plotfig == 1
                    figure;subplot 211;
                    plot(SelPhase1); hold on; plot(Index,SelPhase1(Index),'r*');hold on;
                    subplot 212; plot(x_axis,abs(fftshift(fft(SelPhase1))));
                end
                
            end
            BPM_in_Time = mean(BPM_per_Selphase);
            for ss = 1:SelNum
                SelPhase = [LoadPhase{ChosenIndex(ss)};SelPhase];
            end
            len = length(SelPhase);
            x_axis = [-len/2:len/2-1]*fs(ii)*60/len;
            if plotfig == 1
                    figure;subplot 211;
                    plot(SelPhase); hold on; plot(Index,SelPhase(Index),'r*');hold on;
                    subplot 212; plot(x_axis,abs(fftshift(fft(SelPhase))));
                end
            phaseFre = abs((fft(SelPhase)));
%             [pks, Index, ~, ~] = findpeaks(phaseFre(1:round(len/2)));
            [pks, Index] = sort(phaseFre(1:round(len/2)),'descend');
            LowerIndex = floor(MinBPM/60*len/2/fs(ii)); % 最大呼吸索引
            UpperIndex = ceil(MaxBPM/60*len/2/fs(ii));% 最小呼吸索引
            NormalIndex = round(BPM_in_Time/60*len/2/fs(ii)); % 正常呼吸索引
            MinIndexInterval = ceil(3/60*len/2/fs(ii)); % 正常呼吸索引
            % 离群值及强度联合选取
            IndexLen = 20;
            ChosenIndex = Index(1:IndexLen);
            ChosenIndex = ChosenIndex(ChosenIndex >= LowerIndex);
            ChosenIndex = ChosenIndex(ChosenIndex <= UpperIndex);
            Chosenpks = phaseFre(ChosenIndex);
            Chosenpks = Chosenpks/max(Chosenpks); % 强度归一化
            %         Chosenpks = 1./(1 + exp(Chosenpks)); % 强度归一化
            %         Chosenpks = arctan(Chosenpks*pi/2);
            IndexDiff = abs(ChosenIndex - NormalIndex)/max(abs(ChosenIndex - NormalIndex));
            ChosenWeight = Chosenpks - IndexDiff;
            [~,ChosenSort] = sort(ChosenWeight,'descend');
            ChosenIndex = ChosenIndex(ChosenSort);
            Person = dataPerson(ii);
            for iii = 1:length(ChosenIndex)
                if Person == 0
                    break;
                elseif Person == dataPerson(ii)
                    BPM(ii,Person) = ChosenIndex(iii);
                    Person = Person - 1;
                else
                    if sum(abs(BPM(ii,Person + 1:dataPerson(ii)) - ChosenIndex(iii)) >= MinIndexInterval) == dataPerson(ii) - Person
                        BPM(ii,Person) = ChosenIndex(iii);
                        Person = Person - 1;
                    end
                end
            end
            if Person > 0
                BPM(ii,1:Person) = ChosenIndex(1);
            end
            BPM(ii,:) = BPM(ii,:)*2*fs(ii)/len*60;            
                    
%             while person > 0
%             if Person == dataPerson(ii)
%                 BPM(ii,1) = ChosenIndex(1)*2*fs(ii)/len*60;
%                 Person = Person - 1;
%             else
%                 
%                 BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
%             else
%                 BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
%                 BPM(ii,length(ChosenSort)+1:end) = NormalIndex*2*fs(ii)/len*60;
%             end
%         end

        end
        [ii,BPM_in_Time,BPM(ii,1:dataPerson(ii))]

    end
    
    
    BPM = reshape(BPM',1,[]);
    BPM = BPM(BPM > 0);
    GroundTruth(isnan(GroundTruth)) = 0;
    GroundTruth = reshape(GroundTruth',1,[]);
    GroundTruth = GroundTruth(GroundTruth > 0);
    %
    
    if plotfig == 0
        RMSE_sum = sum((GroundTruth - BPM).^2) + RMSE_sum;
        DataNum = dataNum + DataNum;
        RMSE(SampleNum) = sqrt(1/dataNum*sum((GroundTruth - BPM).^2));
    end
    SampleNum
end

if plotfig == 0
   RMSE
   RMSE_all = sqrt(1/DataNum*RMSE_sum)
end
