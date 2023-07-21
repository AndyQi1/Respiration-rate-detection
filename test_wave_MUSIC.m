%% 导入数据
Y = ['1'];
X = {'001','002','003','004','005','006','007','008','009'};
RMSE = zeros(1,length(X));
RMSE_baseline = zeros(1,length(X));
MinBPM = 5;
MaxBPM = 35;
plotfig = 1;
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
    for ii = 42%:length(N_T)
        
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
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = detrend(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 0); % 去除趋势
                phase(:,(ss-1)*(N_Rx-1) + nn - 1) = hampel(phase(:,(ss-1)*(N_Rx-1) + nn - 1), 10);
        
                Movestd = movstd(phase(:,(ss-1)*(N_Rx-1) + nn - 1),WindowLen);
                Sel = Movestd <= mean(Movestd) + sqrt(var(Movestd));
                StartPoint = sum(Sel(1:round(len/3)) == 0) + 1;
                EndPoint = len - sum(Sel(round(2*len/3):end) == 0);
%                 phase([1:StartPoint,EndPoint:len],(ss-1)*(N_Rx-1) + nn - 1) = 0; % 去除波动较大的值
                LoadPhase{(ss-1)*(N_Rx-1) + nn - 1} = phase([StartPoint:EndPoint],(ss-1)*(N_Rx-1) + nn - 1);
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
        
        ChosenNum = 1;
        MAD = zeros(1, N_Sc*(N_Rx - 1));
        for mm = 1:N_Sc*(N_Rx - 1)
            MAD(mm) = mad(LoadPhase{mm},0,1);
        end
        [~,MADInd] = sort(MAD,'descend');
        ChosenIndex = MADInd(2);
%         Weight = exp(-MAD/mean(MAD));
%         SelPhase = phase.*Weight;
%         SelPhase = reshape(SelPhase,2*len,[]);
%         SelPhase = sum(SelPhase,2);
        SelPhase = LoadPhase{ChosenIndex};
        % 设计滤波器
% fpass = [0.08 0.7]; % 通带频率范围
% 
% % 滤波信号
% SelPhase = bandpass(SelPhase, fpass,fs(ii));
%         N = length(SelPhase);
%         N = N_Sc*N_Rx/2;
%         L = 45;
%         H_ = zeros(L,1);
%         Scalar = 0;
%         H = zeros(N);
%         for tt = 1:N_T(ii)
% %         for hh = 1:N_Sc*N_Rx
%             H = fft(data(tt,1:2:end))'*fft(data(tt,1:2:end)) + H;
%             H = fft(data(tt,2:2:end))'*fft(data(tt,2:2:end)) + H;
% %         end
%         H = 1/(N_T(ii)*N_Sc*N_Rx)*H;
%         end
%         for HH = 1:N + 1-L - 1 % postive smoothing
%         %     eigThreshold = 3e-3;
%         H_ = H_ + H(HH:HH+L-1,HH:HH+L-1);
%         Scalar = Scalar + 1;
%         end
%         F = fliplr(eye(N));                             % transpose matrix
%         H = F*(conj(H))*F;
%         for HH = 1:N + 1-L - 1 % negative smoothing
%         %     eigThreshold = 3e-3;
%         H_ = H_ + H(HH:HH+L-1,HH:HH+L-1);
%         Scalar = Scalar + 1;
%         end
%         H_ = 1/Scalar*H_;
% a = H_;
%     
%     [U,D] = eig(H);
%     %     D = abs(D);
%     D = diag(D)';
%     [D, I] = sort(D);
%     U = fliplr(U(:,I));
%     Res = 1/1000;
%     DelayLen = round(MaxBPM/60/Res);
%     Lp = dataPerson(ii);
%     P_MUSIC = zeros(1,DelayLen);
%     for kk = 1:DelayLen
%         V = exp(-1j*2*pi*[0:length(D)-1]'*delta_f*kk*Res);
%         P_MUSIC(kk) = 1/abs((V'*U(:,Lp + 1:end)*(V'*U(:,Lp + 1:end))'));
%     end
%     %         [v,id] = max(P_MUSIC);
%     %         id*Res/Tc
%     P_MUSIC = 10*log10((P_MUSIC)/max(P_MUSIC));
%     [pks, pksid, w, p] = findpeaks(P_MUSIC);
%     [~,Sortpks] = sort(pks,'descend');
%     pksid(Sortpks(1:Lp))*Res*60
%     BPM(ii,1:dataPerson(ii)) = pksid(Sortpks(1:Lp))*Res*60;
%     BPM(ii,1:dataPerson(ii)) = sort(BPM(ii,1:dataPerson(ii)),'ascend');
%     if plotfig == 1
%         close all;figure;plot(P_MUSIC);
%     end

    end
    
    
    BPM = reshape(BPM',1,[]);
    BPM = BPM(BPM > 0);
    GroundTruth(isnan(GroundTruth)) = 0;
    GroundTruth = reshape(GroundTruth',1,[]);
    GroundTruth = GroundTruth(GroundTruth > 0);
    %
    RMSE(SampleNum) = sqrt(1/dataNum*sum((GroundTruth - BPM).^2));
    SampleNum
end
RMSE
