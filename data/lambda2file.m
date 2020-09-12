%计算数据集各维的lambda，存入文件

filePath =  'Tic-Tac-Toe.xlsx';
[NUM,TXT,RAW] = xlsread(filePath);  %读取Excele表中的数据
[m,n] = size(RAW);

%转换为符号型
RAW = cellfun(@(x){num2str(x)},RAW);

%统计各维符号出现的次数和频率
symbolSta = dataSta( RAW(:,1:n-1));   %dataSta函数统计各统符号频度

result = lambdaD(m, symbolSta, 2);%MSE方法

%存为文件
file_name = ['d:\lambda_', filePath];
xlswrite(file_name, result);