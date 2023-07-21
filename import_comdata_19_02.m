%% 生成比赛提交数据
Y = ['1'];
X = {'001','002','003'};
MinBPM = 12;
MaxBPM = 30;
plotfig = 0;
for SampleNum = 1:length(X)
    %% 导入数据
    name = strcat('Round',Y,'CfgData',X{SampleNum},'.txt');
    Cfgdata = importdata(strcat('CompetitionData1/',name));
    % Cfgdata = importdata('CompetitionData1/Round1CfgData001.txt');
    name = strcat('Round',Y,'Input','Data',X{SampleNum},'.txt');
    Inputdata = importdata(strcat('CompetitionData1/',name));
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
        WindowLen = round(len/(T_Dur(ii)*(MaxBPM + MinBPM)/120)); % 计算滑动窗口
        LoadPhase = cell(1,N_Sc*(N_Rx - 1));
        pr = zeros(1,N_Sc*(N_Rx - 1));
        for ss = 1:N_Sc
            for nn = 2:N_Rx
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = 180*angle(conj(data(:,(ss-1)*N_Rx + 1)).*data(:,(ss-1)*N_Rx + nn))/pi; % 转为角度制
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = detrend(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 1); % 去除趋势
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = hampel(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 10);
        
                Movestd = movstd(phase(:,(ss-1)*(N_Rx-1) + nn - 1),WindowLen);
                Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
                StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
                EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%                 phase([1:StartPoint,EndPoint:len],(ss-1)*(N_Rx-1) + nn - 1) = 0; % 去除波动较大的值
                LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = phase([StartPoint:EndPoint],(ss-1)*(N_Rx-1) + nn - 1);
%                 fpass = [0.1 0.6]; % 通带频率范围
% % % % % 滤波信号
% LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = bandpass(LoadPhase{(ss-1)*(N_Rx-1) + nn - 1}, fpass,fs(ii));
%                 % 基于重现图选取子载波
%                 RPLen = 100;
%                 RP = zeros(RPLen);
%                 epsilon = 10;
%                 for rr = 1:RPLen
%                     for pp = 1:RPLen
%                         if norm(LoadPhase{(ss-1)*(N_Rx-1) + nn - 1}(round((EndPoint - StartPoint)/2) - RPLen/2 + rr)...
%                                 - LoadPhase{(ss-1)*(N_Rx-1) + nn - 1}(round((EndPoint - StartPoint)/2) - RPLen/2 + pp)) < epsilon
%                             RP(rr,pp) = 1;
%                         end
%                     end
%                 end
%                 RP = imrotate(RP,45);
%                 sigval = svd(RP);
%                 pr((ss-1)*(N_Rx-1) + nn - 1) = sigval(1)/sum(sigval);
            end
            
        end
        
        % 子载波筛选
%         [~,ChosenIndex] = max(pr);
        ChosenStart = 1;
        ChosenNum = 10;
        MAD = zeros(1, N_Sc*(N_Rx - 1));
        for mm = 1:N_Sc*(N_Rx - 1)
            MAD(mm) = mad(LoadPhase{mm},0);
        end
        fpass = [0.1 0.6]; % 通带频率范围
        % % % % 滤波信号
        
        [~,MADInd] = sort(MAD,'descend');
        pkssum = zeros(1,ChosenNum);
        for cc = ChosenStart:ChosenStart + ChosenNum - 1
            LoadPhase{MADInd(cc)} = bandpass(LoadPhase{MADInd(cc)}, fpass, fs(ii));
            Phasecorr = xcorr(LoadPhase{MADInd(cc)});
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
        SelPhase = [];
        for ss = 1:2
                SelPhase = [LoadPhase{ChosenIndex(ss)};SelPhase];
            end
%         if dataPerson(ii) == 1
%             SelPhase = LoadPhase{ChosenIndex(1)};
%         elseif dataPerson(ii) == 2
%             for ss = 1:3
%                 SelPhase = [LoadPhase{ChosenIndex(ss)};SelPhase];
%             end
%         else
%             for ss = 1:3
%                 SelPhase = [LoadPhase{ChosenIndex(ss)};SelPhase];
%             end
%         end
        
        len = length(SelPhase);
%         SelPhase = phase(:,MADInd(7));
%         SelPhase = hampel(SelPhase);
        SelPhase = reshape(SelPhase,[],1);
        x_axis = [-len/2:len/2-1]*fs(ii)*60/len;
        if plotfig == 1
            close all;figure;subplot 211; plot(SelPhase);
            subplot 212; plot(x_axis,abs(fftshift(fft(SelPhase))));
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
            figure;subplot 211; plot(SelPhase);
            subplot 212; plot(x_axis,abs(fftshift(fft(SelPhase))));
        end
        %         hpf = designfilt('highpassfir', 'FilterOrder', 100, 'CutoffFrequency', 1, 'SampleRate', fs(ii));
        phaseFre = abs((fft(SelPhase)));
%         close all
%         semilogy(phaseFre(1:round(len/2)));
%         len = length(phaseFre);
%         [~,Index] = sort(phaseFre(1:round(len/2)),'descend');
        [pks, Index, w, p] = findpeaks(phaseFre(1:round(len/2)));
        LowerIndex = ceil(MinBPM/60*len/2/fs(ii)); % 最大呼吸索引
        UpperIndex = floor(MaxBPM/60*len/2/fs(ii));% 最小呼吸索引
        NormalIndex = round(20/60*len/2/fs(ii)); % 正常呼吸索引
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
%         Chosenpks = Chosenpks/max(Chosenpks); % 强度归一化
        Chosenpks = 1./(1 + exp(Chosenpks)); % 强度归一化
        IndexDiff = abs(ChosenIndex - NormalIndex)/max(abs(ChosenIndex - NormalIndex));
        ChosenWeight = Chosenpks - IndexDiff;
        [~,ChosenSort] = sort(ChosenWeight,'descend');
        if length(ChosenSort) >= dataPerson(ii)
            BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
        else
            BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
            BPM(ii,length(ChosenSort):end) = NormalIndex*2*fs(ii)/len*60;
        end
        [ii,BPM(ii,1:dataPerson(ii))]
%         BPM(ii,1:dataPerson(ii)) = FindBPM*2*fs(ii)/len*60;
% if dataPerson(ii) == 1
%         BPM(ii,1:dataPerson(ii))  = 20 + randn(1,dataPerson(ii)); % 随机产生数据
% elseif dataPerson(ii) == 2
%     BPM(ii,1) = 10 + randn;
%     BPM(ii,2) = 20 + randn;
% else
%     BPM(ii,1) = 10 + randn;
%     BPM(ii,2) = 15 + randn;
%     BPM(ii,3) = 20 + randn;
% end
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
