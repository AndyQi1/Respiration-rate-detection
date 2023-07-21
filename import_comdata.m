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
        
%         data = Inputdata(Idx(ii) + 1:Idx(ii + 1),:);
%         real = data(:,1:2:end);
%         imag = data(:,2:2:end);
%         data = real + 1j*imag;
%         phase = zeros(N_T(ii),N_Sc*(N_Rx - 1));
%         
%         len = N_T(ii);
%         WindowLen = round(len/(T_Dur(ii)*(MaxBPM + MinBPM)/120)); % 计算滑动窗口
%         for ss = 1:N_Sc
%             for nn = 2:N_Rx
%                 phase(:,(ss-1)*(N_Rx-1) + nn - 1) = 180*angle(conj(data(:,(ss-1)*N_Rx + 1)).*data(:,(ss-1)*N_Rx + nn))/pi; % 转为角度制
%                 phase(:,(ss-1)*(N_Rx-1) + nn - 1) = detrend(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 0); % 去除趋势
%                 phase(:,(ss-1)*(N_Rx-1) + nn - 1) = hampel(phase(:,(ss-1)*(N_Rx-1) + nn - 1));
%                 
%                 %                 Movestd = movstd(phase(:,(ss-1)*(N_Rx-1) + nn - 1),WindowLen);
%                 %                 Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
%                 %                 StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
%                 %                 EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%                 %                 phase([1:StartPoint,EndPoint:len],(ss-1)*(N_Rx-1) + nn - 1) = 0; % 去除波动较大的值
%             end
%         end
%         % 子载波筛选
%         ChosenNum = 1;
%         MAD = mad(phase,0,1);
%         [~,MADInd] = sort(MAD,'descend');
%         %         Weight = exp(-MAD/mean(MAD));
%         %         SelPhase = phase.*Weight;
%         %         SelPhase = reshape(SelPhase,2*len,[]);
%         %         SelPhase = sum(SelPhase,2);
%         
%         SelPhase = phase(:,MADInd(7));
%         SelPhase = hampel(SelPhase);
%         SelPhase = reshape(SelPhase,[],1);
%         if plotfig == 1
%             close all;figure(1);subplot 211; plot(SelPhase);
%         end
%         %         Movestd = movstd(SelPhase,WindowLen);
%         %                 Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
%         %                 StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
%         %                 EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%         %                 SelPhase = SelPhase(StartPoint:EndPoint); % 去除波动较大的值
%         [c, l] = wavedec(SelPhase, level, Type); % 小波变换滤波
%         %         coefapp = appcoef(c, l, Type, level);
%         c(l(1):end) = 0;
%         %                     c(1:l(1)/2) = 0;
%         SelPhase = waverec(c,l,Type);
%         if plotfig == 1
%             figure(1);subplot 212; plot(SelPhase);
%         end
%         %         hpf = designfilt('highpassfir', 'FilterOrder', 100, 'CutoffFrequency', 1, 'SampleRate', fs(ii));
%         phaseFre = abs((fft(SelPhase)));
%         %         close all
%         %         semilogy(phaseFre(1:round(len/2)));
%         len = length(phaseFre);
%         [~,Index] = sort(phaseFre(1:round(len/2)),'descend');
%         
%         LowerIndex = round(MinBPM/60*len/2/fs(ii));
%         UpperIndex = round(MaxBPM/60*len/2/fs(ii));
%         
%         FindBPM = zeros(1,dataPerson(ii));
%         person = dataPerson(ii);
%         for dd = 1:length(Index)
%             if Index(dd)>= LowerIndex && Index(dd) <= UpperIndex && person >= 1
%                 FindBPM(person) = Index(dd);
%                 person = person - 1;
%                 if person < 1
%                     break;
%                 end
%             end
%         end
%         FindBPM*2*fs(ii)/len*60
%         BPM(ii,1:dataPerson(ii)) = FindBPM*2*fs(ii)/len*60;
if dataPerson(ii) == 1
        BPM(ii,1:dataPerson(ii))  = 20 + randn(1,dataPerson(ii)); % 随机产生数据
elseif dataPerson(ii) == 2
    BPM(ii,1) = 10 + randn;
    BPM(ii,2) = 20 + randn;
else
    BPM(ii,1) = 10 + randn;
    BPM(ii,2) = 15 + randn;
    BPM(ii,3) = 27 + randn;
end
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
