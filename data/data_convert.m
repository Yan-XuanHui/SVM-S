%����ֵ������ת��Ϊ��������
clear;
[NUM,TXT,RAW]=xlsread('N102400A8.xlsx');
%[NUM,TXT,RAW]=csvread('N100A8.csv');
[n,m] = size(RAW);

%ȡ��1������(��ֵ�ͣ�
for col=1:m  %��
    if  ( strcmp(class(cell2mat(RAW(1,col))),'double') )
        aa = cell2mat(RAW(:,col));  %��һ������
        mina = min(aa);
        maxa = max(aa);
        bins = 10;  %����Ϊ10��
        char = ['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V'];
        size_aa = size(aa,1);
        for i = 1:size_aa  %��
            k = (aa(i) - mina) / ((maxa - mina)/bins)+1  %�����ڵڼ���bin
            RAW{i,col} = char( int16(k) ,1);
        end; 
    end;
    
end;


s = xlswrite('N102400A8_C.xlsx', RAW); 