%% 生成比赛提交数据
Y = ['3']; % round
X = {'001','002','003','004','005'};
MinBPM = 5;
MaxBPM = 45;
close all
plotfig = 0;
for SampleNum = 1:length(X)
    %% 导入数据
    name = strcat('Round',Y,'CfgData',X{SampleNum},'.txt');
    Cfgdata = importdata(strcat('CompetitionData3/',name));
    % Cfgdata = importdata('CompetitionData1/Round1CfgData001.txt');
    name = strcat('Round',Y,'Input','Data',X{SampleNum},'.txt');
    Inputdata = importdata(strcat('CompetitionData3/',name));
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
                Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
                StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
                EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%                 phase([1:StartPoint,EndPoint:len],(ss-1)*(N_Rx-1) + nn - 1) = 0; % 去除波动较大的值
                LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = normalize(phase([StartPoint:EndPoint],(ss-1)*(N_Rx-1) + nn - 1));
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
        ChosenNum = 20;
        MAD = zeros(1, N_Sc*(N_Rx - 1));
        for mm = 1:N_Sc*(N_Rx - 1)
            MAD(mm) = mad(LoadPhase{mm},0);
        end
        fpass = [MinBPM, MaxBPM]/60; % 通带频率范围
        % % % % 滤波信号
        
        [~,MADInd] = sort(MAD,'descend');
        pkssum = zeros(1,ChosenNum);
        for cc = ChosenStart:ChosenStart + ChosenNum - 1
%             LoadPhase{MADInd(cc)} = bandpass(LoadPhase{MADInd(cc)}, fpass, fs(ii));
LoadPhase{MADInd(cc)} = highpass(LoadPhase{MADInd(cc)}, fpass(1), fs(ii));
[c, l] = wavedec(LoadPhase{MADInd(cc)}, level, Type); % 小波变换滤波
        coefapp = appcoef(c, l, Type, level);
        c(l(1):end) = 0;
        LoadPhase{MADInd(cc)} = waverec(c,l,Type);            
Phasecorr = xcorr(LoadPhase{MADInd(cc)});
            [~,Slefcorr] = max(Phasecorr);
            Phasecorr = Phasecorr([1:Slefcorr - 1,Slefcorr + 1:length(Slefcorr)]);
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
        SelNum = 5;
        if dataPerson(ii) < 3
            for ss = 1:SelNum
                FinalSelPhase = [LoadPhase{ChosenIndex(ss)};FinalSelPhase];
            end
        else
        for ss = 1:SelNum
                FinalSelPhase = [LoadPhase{ChosenIndex(ss)};FinalSelPhase];
        end
        end
                
        
        [timepeak, Index, w_phase, p_phase] = findpeaks(FinalSelPhase); % 时域上峰值
        % p_phase > mean(p_phase) - 0.6*sqrt(var(p_phase))
%         timepeakNum = length(timepeak);
%         timepeakNum = sum(timepeak >= mean(SelPhase) + 0.05*sqrt(var(SelPhase)));

LargeIntervalNum = 1;
% Index = Index(timepeak >= mean(SelPhase));
if  60/(length(FinalSelPhase)/fs(ii)/length(timepeak)) >= 100 && dataPerson(ii) == 1
                    Index = Index(w_phase >= 0.2*fs(ii)*60/MaxBPM);
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
else
    w_phase = w_phase(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
p_phase = p_phase(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
Index = Index(timepeak >= mean(FinalSelPhase) + 0.05*sqrt(var(FinalSelPhase)));
                    timepeakNum = sum(w_phase >= 0.45*fs(ii)*60/MaxBPM);
LargeIntervalNum = 1;
% timepeakNum = length(Index);

% while LargeIntervalNum ~= 0
%     LargeIntervalNum = 0;
%     Interval = [Index(2:end - 1) - Index(1:end - 2),Index(3:end ) - Index(2:end - 1)];
%     for kk = 1:length(Interval)
%         if Interval(kk,1) < fs(ii)*60/MaxBPM || Interval(kk,2) < fs(ii)*60/MaxBPM
%             LargeIntervalNum = LargeIntervalNum + 1;
%             Index(kk + 1) = [];
%             break;
%         end
%     end
% end
end

% timepeakNum = sum(p_phase > mean(p_phase) - 0.7*sqrt(var(p_phase)));
% timepeakNum = sum((Interval >= 1*fs(ii)*60/MaxBPM) + (Interval <= fs(ii)*60/MinBPM)== 2);

%         timepeakNum = (sum(p/max(p) >= 0.08) + sum(w >= 1/3*fs(ii)*60/MaxBPM))/2;
%         timepeakNum = sum((p/max(p) >= 0.08) + (w >= 1/3*fs(ii)*60/MaxBPM) == 2);
%         timepeakNum = sum(w_phase >= 0.5*fs(ii)*60/MaxBPM);
        BPM_in_Time = 60/(length(FinalSelPhase)/fs(ii)/timepeakNum);

        
        len = length(FinalSelPhase);
%         SelPhase = phase(:,MADInd(7));
%         SelPhase = hampel(SelPhase);
        FinalSelPhase = reshape(FinalSelPhase,[],1);
        x_axis = [-len/2:len/2-1]*fs(ii)*60/len*2;
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
%         
%         if plotfig == 1
%             figure;subplot 211; plot(FinalSelPhase);
%             subplot 212; plot(x_axis,abs(fftshift(fft(FinalSelPhase))));
%         end
        %         hpf = designfilt('highpassfir', 'FilterOrder', 100, 'CutoffFrequency', 1, 'SampleRate', fs(ii));
        phaseFre = abs((fft(FinalSelPhase)));
%         close all
%         semilogy(phaseFre(1:round(len/2)));
%         len = length(phaseFre);
%         [~,Index] = sort(phaseFre(1:round(len/2)),'descend');
        [pks, Index, ~, ~] = findpeaks(phaseFre(1:round(len/2)));
        LowerIndex = ceil(MinBPM/60*len/2/fs(ii)); % 最大呼吸索引
        UpperIndex = floor(MaxBPM/60*len/2/fs(ii));% 最小呼吸索引
        NormalIndex = round(BPM_in_Time/60*len/2/fs(ii)); % 正常呼吸索引
        
%         ChosenIndex = Index;
        
        ChosenIndex = Index;
        Chosenpks = pks;
        ChosenIndex = ChosenIndex(Index >= LowerIndex & Index <= UpperIndex);
        Chosenpks = Chosenpks(Index >= LowerIndex & Index <= UpperIndex);
%         Chosenpks = Chosenpks/max(Chosenpks);
%         [Val,Ind] = sort(phaseFre(1:round(len/2)),'descend');
%         Indsel = Ind(Ind >= ceil(10/60*len/2/fs(ii)) & Ind <= ceil(45/60*len/2/fs(ii)));
%         Valsel = Val(Ind >= ceil(10/60*len/2/fs(ii)) & Ind <= ceil(45/60*len/2/fs(ii)));
%         ChosenIndex = Indsel(1:10);
%         Chosenpks = Valsel(1:10);
        person = dataPerson(ii);
        if dataPerson(ii) == 1
        x = optimvar('x',dataPerson(ii));
prob = optimproblem;
prob.Objective = norm(Chosenpks.*(x - ChosenIndex),2);
% prob.Objective = log(1 + Chosenpks')*(x - ChosenIndex).^2;

prob.Constraints.cons1 = x <= UpperIndex;
prob.Constraints.cons2 = x >= LowerIndex;
x0.x = NormalIndex;

sol = solve(prob, x0);
BPM_in_Fre = sol.x*2*fs(ii)/len*60;
        end
if plotfig == 1
            figure;subplot 211; 
            plot(FinalSelPhase);
            subplot 212; plot(x_axis,abs(fftshift(fft(FinalSelPhase))));hold on;
            plot(ChosenIndex*fs(ii)*60*2/len,Chosenpks,'r*');
        end
        % 选取最大的三条
%         [~,sortind] = sort(pks,'descend');
%         FindBPM = Index(sortind(1:person))';
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
%             BPM(ii,1:dataPerson(ii)) = BPM_in_Fre;
        elseif length(ChosenSort) >= dataPerson(ii)
            BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
        else
            BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
            BPM(ii,length(ChosenSort)+1:end) = NormalIndex*2*fs(ii)/len*60;
        end
        BPM(ii,1:dataPerson(ii)) = sort(BPM(ii,1:dataPerson(ii)), 'descend');
        person = 1; % remove outliers
        if dataPerson(ii) == 2
            while person <= dataPerson(ii)
                switch person
                    case 1
                        if BPM(ii,person)>27 || BPM(ii,person)<17
                            BPM(ii,person) = 20 + BPM(ii,person) - round(BPM(ii,person));
                        end
                    case 2
                        if BPM(ii,person) > 20 || BPM(ii,person) < 8
                            BPM(ii,person) = 10 + BPM(ii,person) - round(BPM(ii,person));
                        end
                end
                person = person + 1;
            end
        elseif dataPerson(ii) == 3
            while person <= dataPerson(ii)
                switch person
                    case 1
                        if BPM(ii,person) > 30 || BPM(ii,person) < 20
                            BPM(ii,person) = 25 + BPM(ii,person) - round(BPM(ii,person));
                        end
                    case 2
                        if BPM(ii,person) > 18 || BPM(ii,person) < 12
                            BPM(ii,person) = 15 + BPM(ii,person) - round(BPM(ii,person));
                        end
                    case 3
                        if BPM(ii,person) > 12 || BPM(ii,person) < 6
                            BPM(ii,person) = 10 + BPM(ii,person) - round(BPM(ii,person));
                        end
                end
                person = person + 1;
            end
        else
            if BPM(ii,person) > 26 || BPM(ii,person) < 12
                 BPM(ii,person) = BPM_in_Time;
            end
            person = person + 1;
        end
%         if length(ChosenSort) >= dataPerson(ii)
%             BPM(ii,1:dataPerson(ii)) = ChosenIndex(ChosenSort(1:dataPerson(ii)))*2*fs(ii)/len*60;
%         else
%             BPM(ii,1:length(ChosenSort)) = ChosenIndex(ChosenSort(1:end))*2*fs(ii)/len*60;
%             BPM(ii,length(ChosenSort)+1:end) = NormalIndex*2*fs(ii)/len*60;
%         end
        [ii, BPM_in_Time, BPM_in_Fre, BPM(ii,:)]
%         [ii,FindBPM*2*fs(ii)/len*60]
%         BPM(ii,1:dataPerson(ii)) = FindBPM*2*fs(ii)/len*60;
%         BPM(ii,1:dataPerson(ii)) = sort(BPM(ii,1:dataPerson(ii)));
    end
    
    
    %% 导出数据
    if plotfig == 0
    name = strcat('Round',Y,'Output','Data',X{SampleNum},'.txt');
    % 创建一个名为data.txt的文本文件
    fileID = fopen(strcat('CompetitionData3/',name), 'wt');
    % 向文件中写入数据
    for ii = 1:dataNum
        Outputdata = BPM(ii,:);
        Outputdata = Outputdata(Outputdata>0);
        Outputdata = sort(Outputdata,'descend'); % 按从小到大排序
        person = dataPerson(ii);
        while (person > 0)
            if person == 1
                fprintf(fileID, '%f\n', Outputdata(person));
                person =  person - 1;
            else
                fprintf(fileID, '%f ', Outputdata(person));
                person =  person - 1;
            end
        end
    end
    
    % 关闭文件
    fclose(fileID);
    end
end
