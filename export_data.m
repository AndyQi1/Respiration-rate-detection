%% 输出数据 (测试)
Y = ['1'];
X = {'001','002','003'};
Num = [67,50,58];

for SampleNum = 1:length(X)
    name = strcat('Round',Y,'CfgData',X{SampleNum},'.txt');
    Cfgdata = importdata(strcat('CompetitionData1/',name));
    dataNum = Cfgdata(1); % 数据数
    dataPerson = Cfgdata(2:dataNum + 1); % 每条数据的人数
    name = strcat('Round',Y,'Output','Data',X{SampleNum},'.txt');
    % 创建一个名为data.txt的文本文件
    fileID = fopen(strcat('CompetitionData1/',name), 'w');

    % 向文件中写入数据
    for ii = 1:dataNum
        data = rand(1,dataPerson(ii))*(45 - 5) + 5;
        data = sort(data,'descend');
        while (dataPerson(ii) > 0)
            if dataPerson(ii) == 1
                fprintf(fileID, '%f\n', data(dataPerson(ii)));
                dataPerson(ii) =  dataPerson(ii) - 1;
            else
                fprintf(fileID, '%f ', data(dataPerson(ii)));
                dataPerson(ii) =  dataPerson(ii) - 1;
            end
        end
    end

    % 关闭文件
    fclose(fileID);
end
