%�������ݼ���ά��lambda�������ļ�

filePath =  'Tic-Tac-Toe.xlsx';
[NUM,TXT,RAW] = xlsread(filePath);  %��ȡExcele���е�����
[m,n] = size(RAW);

%ת��Ϊ������
RAW = cellfun(@(x){num2str(x)},RAW);

%ͳ�Ƹ�ά���ų��ֵĴ�����Ƶ��
symbolSta = dataSta( RAW(:,1:n-1));   %dataSta����ͳ�Ƹ�ͳ����Ƶ��

result = lambdaD(m, symbolSta, 2);%MSE����

%��Ϊ�ļ�
file_name = ['d:\lambda_', filePath];
xlswrite(file_name, result);