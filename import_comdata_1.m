%% 生成比赛提交数据
Y = ['2']; % round
X = {'001','002','003','004'}; 
MinBPM = 5;
MaxBPM = 45;
plotfig = 0;
for SampleNum = 1:length(X)
    %% 导入数据
    name = strcat('Round',Y,'CfgData',X{SampleNum},'.txt');
    Cfgdata = importdata(strcat('CompetitionData2/',name));
    % Cfgdata = importdata('CompetitionData1/Round1CfgData001.txt');
    name = strcat('Round',Y,'Input','Data',X{SampleNum},'.txt');
    Inputdata = importdata(strcat('CompetitionData2/',name));
    % Breathwave = importdata('TestData/TestBreathWave001.txt');
    
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
    % 数据每列先遍历接收天线，再遍历子载波， 每行代表时序上的信息
    Idx = cumsum([0;N_T]);
    BPM = zeros(dataNum,3);
    %% 数据处理
   for ii = 1:length(N_T)
        
        data = Inputdata(Idx(ii) + 1:Idx(ii + 1),:);
        real = data(:,1:2:end);
        imag = data(:,2:2:end);
        data = real + 1j*imag;
        phase = zeros(N_T(ii),N_Sc*(N_Rx - 1));
        
        len = N_T(ii);
        Precyclix = round(3/(T_Dur(ii)/len)); % 保护区间（5s）
        WindowLen = round(len/(T_Dur(ii)*(MaxBPM + MinBPM)/120)); % 计算滑动窗口
        LoadPhase = cell(1,N_Sc*(N_Rx - 1));
        pr = zeros(1,N_Sc*(N_Rx - 1));
        for ss = 1:N_Sc
            for nn = 2:N_Rx
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = 180*angle(conj(data(:,(ss-1)*N_Rx + 1)).*data(:,(ss-1)*N_Rx + nn))/pi; % 转为角度制
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = detrend(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 1); % 去除趋势
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = hampel(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 10);
        
                Movestd = movstd(phase(:,(ss-1)*(N_Rx-1) + nn - 1),WindowLen);
%                 Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
%                 StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
%                 EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
                
%                 phase([1:StartPoint,EndPoint:len],(ss-1)*(N_Rx-1) + nn - 1) = 0; % 去除波动较大的值
%                 StartPoint = 1;EndPoint = len;Precyclix = 0;
%                 LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = phase([StartPoint + Precyclix:EndPoint - Precyclix]...
%                     ,(ss-1)*(N_Rx-1) + nn - 1);
%                 
%                 LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = normalize(phase([StartPoint + Precyclix:EndPoint - Precyclix]...
%                     ,(ss-1)*(N_Rx-1) + nn - 1));
                if plotfig == 1
                    figure(1);plot(phase(:,(ss-1)*(N_Rx-1) + nn - 1));hold on;
                end
%                 
            end
            
        end
        SelMagNum = 10;
        phase = phase(Precyclix:end - Precyclix,:); % 保护间隔
        [~,SortMag] = sort(max(abs(phase)),'descend'); % 根据幅度排序
        SelPhase = phase(:,SortMag(1:SelMagNum));
        SelPhase = normalize(SelPhase);
        % 子载波筛选
%         [~,ChosenIndex] = max(pr);
        ChosenStart = 1;
        ChosenNum = 5;
        MAD = zeros(1, SelMagNum);
        for mm = 1:SelMagNum
            MAD(mm) = mad(SelPhase(:,mm),0);
        end
        fpass = [MinBPM, MaxBPM]/60; % 通带频率范围
        % % % % 滤波信号
        
        [~,MADInd] = sort(MAD,'descend');
        pkssum = zeros(1,ChosenNum);
        for cc = ChosenStart:ChosenStart + ChosenNum - 1
            SelPhase(:,MADInd(cc)) = bandpass(SelPhase(:,MADInd(cc)), fpass, fs(ii));
            Phasecorr = xcorr(SelPhase(:,MADInd(cc)));
            [~,Selfcorr] = max(Phasecorr);
            Phasecorr = Phasecorr([1:Selfcorr - 1,Selfcorr + 1:length(Selfcorr)]);
            Phasecorr = Phasecorr/max(Phasecorr);
            [pks, ~, ~, ~] = findpeaks(Phasecorr);
            pks = sort(pks,'descend');
            pkssum(cc) = sum(pks(1:min(length(pks),10)));
        end
        [~,pkssumIndex] = sort(pkssum,'descend');
        
%         ChosenIndex = MADInd(1);
%         Weight = exp(-MAD/mean(MAD));
%         SelPhase = phase.*Weight;
%         SelPhase = reshape(SelPhase,2*len,[]);
%         SelPhase = sum(SelPhase,2);
ChosenIndex = MADInd(pkssumIndex(1:ChosenNum)); % 选取三条
        FinalSelPhase = [];
        SelNum = 2;
        if dataPerson(ii) < 3
            for ss = 1:SelNum
                FinalSelPhase = [SelPhase(:,ChosenIndex(ss));FinalSelPhase];
            end
        else
        for ss = 1:SelNum
                FinalSelPhase = [SelPhase(:,ChosenIndex(ss));FinalSelPhase];
        end
        end
        [timepeak, Index, w_phase, p_phase] = findpeaks(FinalSelPhase); % 时域上峰值
        % p_phase > mean(p_phase) - 0.6*sqrt(var(p_phase))
%         timepeakNum = length(timepeak);
%         timepeakNum = sum(timepeak >= mean(SelPhase) + 0.05*sqrt(var(SelPhase)));
w_phase = w_phase(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
p_phase = p_phase(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
Index = Index(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
% Index = Index(timepeak >= mean(SelPhase));
LargeIntervalNum = 1;

while LargeIntervalNum ~= 0
    LargeIntervalNum = 0;
    Interval = [Index(2:end - 1) - Index(1:end - 2),Index(3:end ) - Index(2:end - 1)];
    for kk = 1:length(Interval)
        if Interval(kk,1) < fs(ii)*60/MaxBPM || Interval(kk,2) < fs(ii)*60/MaxBPM
            LargeIntervalNum = LargeIntervalNum + 1;
            Index(kk + 1) = [];
            break;
        end
    end
end
timepeakNum = length(Index);
% timepeakNum = sum(p_phase > mean(p_phase) - 0.7*sqrt(var(p_phase)));
% timepeakNum = sum((Interval >= 1*fs(ii)*60/MaxBPM) + (Interval <= fs(ii)*60/MinBPM)== 2);

%         timepeakNum = (sum(p/max(p) >= 0.08) + sum(w >= 1/3*fs(ii)*60/MaxBPM))/2;
%         timepeakNum = sum((p/max(p) >= 0.08) + (w >= 1/3*fs(ii)*60/MaxBPM) == 2);
%         timepeakNum = sum(w >= 0.5*fs(ii)*60/MaxBPM);
        BPM_in_Time = 60/(length(FinalSelPhase)/fs(ii)/timepeakNum);

        
        len = length(FinalSelPhase);
%         SelPhase = phase(:,MADInd(7));
%         SelPhase = hampel(SelPhase);
        FinalSelPhase = reshape(FinalSelPhase,[],1);
        x_axis = [-len/2:len/2-1]*fs(ii)*60/len;
        if plotfig == 1
            figure;subplot 211; plot(Index,FinalSelPhase(Index),'r*');hold on;
            plot(FinalSelPhase);
            subplot 212; plot(x_axis,abs(fftshift(fft(FinalSelPhase))));
        end
%         Movestd = movstd(SelPhase,WindowLen);
%                 Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
%                 StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
%                 EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%                 SelPhase = SelPhase(StartPoint:EndPoint); % 去除波动较大的值
% fpass = [0.08 0.65]; % 通带频率范围
% fpass = [0.1 0.6]; % 通带频率范围
% % % % 滤波信号
% SelPhase = bandpass(SelPhase, fpass,fs(ii));

%         [c, l] = wavedec(SelPhase, level, Type); % 小波变换滤波
%         coefapp = appcoef(c, l, Type, level);
%         c(l(1):end) = 0;
%         SelPhase = waverec(c,l,Type);
%         SelPhase = sgolayfilt(SelPhase, 5, 101);
        
        if plotfig == 1
            figure;subplot 211; plot(FinalSelPhase);
            subplot 212; plot(x_axis,abs(fftshift(fft(FinalSelPhase))));
        end
        %         hpf = designfilt('highpassfir', 'FilterOrder', 100, 'CutoffFrequency', 1, 'SampleRate', fs(ii));
        phaseFre = abs((fft(FinalSelPhase)));
%         close all
%         semilogy(phaseFre(1:round(len/2)));
%         len = length(phaseFre);
%         [~,Index] = sort(phaseFre(1:round(len/2)),'descend');
        [pks, Index, ~, ~] = findpeaks(phaseFre(1:round(len/2)));
        LowerIndex = floor(MinBPM/60*len/2/fs(ii)); % 最大呼吸索引
        UpperIndex = ceil(MaxBPM/60*len/2/fs(ii));% 最小呼吸索引
        NormalIndex = round(BPM_in_Time/60*len/2/fs(ii)); % 正常呼吸索引
        FindBPM = zeros(1,dataPerson(ii));
        person = dataPerson(ii);
        % 选取最大的三条
%         [~,sortind] = sort(pks,'descend');
%         FindBPM = Index(sortind(1:person))';
        % 扫描选取
        for dd = 1:length(Index)
            if Index(dd)>= LowerIndex && Index(dd) <= UpperIndex && person >= 1
                FindBPM(person) = Index(dd);
                person = person - 1;
                if person < 1
                    break;
                end
            end
        end
        if person >= 1
            FindBPM(1:person) = Index(1:person);
        end
        % 离群值及强度联合选取
        ChosenIndex = Index;
        ChosenIndex = ChosenIndex(ChosenIndex >= LowerIndex);
        ChosenIndex = ChosenIndex(ChosenIndex <= UpperIndex);
        Chosenpks = phaseFre(ChosenIndex);
        Chosenpks = Chosenpks/max(Chosenpks); % 强度归一化
%         Chosenpks = 1./(1 + exp(Chosenpks)); % 强度归一化
%         Chosenpks = arctan(Chosenpks*pi/2);
        IndexDiff = abs(ChosenIndex - NormalIndex)/max(abs(ChosenIndex - NormalIndex));
        ChosenWeight = Chosenpks - IndexDiff;
        [~,ChosenSort] = sort(ChosenWeight,'descend');
        if dataPerson(ii) == 1
            BPM(ii,1:dataPerson(ii)) = (BPM_in_Time + ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60)/2;
            BPM(ii,1:dataPerson(ii)) = BPM_in_Time;
        elseif length(ChosenSort) >= dataPerson(ii)
            BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
        else
            BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
            BPM(ii,length(ChosenSort)+1:end) = NormalIndex*2*fs(ii)/len*60;
        end
%         if length(ChosenSort) >= dataPerson(ii)
%             BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
%         else
%             BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
%             BPM(ii,length(ChosenSort)+1:end) = NormalIndex*2*fs(ii)/len*60;
%         end
        [ii,BPM_in_Time,BPM(ii,1:dataPerson(ii))]
%         [ii,FindBPM*2*fs(ii)/len*60]
%         BPM(ii,1:dataPerson(ii)) = FindBPM*2*fs(ii)/len*60;
%         BPM(ii,1:dataPerson(ii)) = sort(BPM(ii,1:dataPerson(ii)));
    end
    
    
    %% 导出数据
    name = strcat('Round',Y,'Output','Data',X{SampleNum},'.txt');
    % 创建一个名为data.txt的文本文件
    fileID = fopen(strcat('CompetitionData1/',name), 'w');
    % 向文件中写入数据
    for ii = 1:dataNum
        Outputdata = BPM(ii,:);
        Outputdata = Outputdata(Outputdata>0);
        Outputdata = sort(Outputdata,'descend'); % 按从小到大排序
        while (dataPerson(ii) > 0)
            if dataPerson(ii) == 1
                fprintf(fileID, '%f\n', Outputdata(dataPerson(ii)));
                dataPerson(ii) =  dataPerson(ii) - 1;
            else
                fprintf(fileID, '%f ', Outputdata(dataPerson(ii)));
                dataPerson(ii) =  dataPerson(ii) - 1;
            end
        end
    end
    
    % 关闭文件
    fclose(fileID);
end
