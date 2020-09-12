function [result]=SVM_Categoricla_Fun(dataSet, K, KO1, KO2)  
% INPUT
% dataSet: Dataset's name，
% K: Kernel (gaussian or poly)，
% KO1: Parameter of SVM (degree or gamma)
% KO2: Parameter C of Gaussian Kernel
%%

filePath =  dataSet;
fprintf('Kernel:%s, KernelOption:%f\n', K, KO1);

[NUM,TXT,RAW] = xlsread(filePath);  %read data from Excel file
%[NUM,TXT,RAW] = csvread('GermanCredit.xlsx')
[m,n] = size(RAW);

%size of train set
m1 = round(m*0.7);
%size of test set
m2 = m - m1;

%convert numerical data to symbolic data
RAW = cellfun(@(x){num2str(x)},RAW);     

trainSet =  RAW(randperm(m, m1),:) ;  % randomly generate the train set
testSet = RAW(randperm(m, m2),:) ;

%Count the data frequency of each attribute
symbolSta = dataSta( RAW(:,1:n-1) );   %dataSt()a: function of data frequency statistic

%compute lambda
my_lambda = lambdaD(m, symbolSta, 2); % computer bandwidth with MSE method

%create train data
xapp = trainSet(:,1:n-1);
yapp_RAW = trainSet(:,n);
yapp_sta = tabulate(yapp_RAW);
yapp = zeros(m1,1);

%convert label to -1 or 1
for i=1:m1
    if strcmp(yapp_sta{1}, yapp_RAW(i))
        yapp(i) = 1;
    else
        yapp(i) = -1;
    end;
end;


%%
%training 
kernel = lower(K);   % lower-case
kerneloption = KO1;
lambda = 1e-7;  
C = KO2;         %bound on lagrangian multiplier

if strcmp(kernel,'poly')   
    [dotProduct] = dotProductMatrix(xapp, symbolSta, my_lambda);  %compute dot product and store in the matrix
else
    dotProduct = ones(1,1);    % dot product is computed in svmkernel_categorical.m
end;

[xsup,w,w0,pos,tps,alpha] = svmclass_categorical(xapp,yapp, my_lambda, C,lambda,kernel,kerneloption,1, symbolSta, dotProduct); 


%%
%testing
xtest = testSet(:,1:n-1);
ytest_RAW = testSet(:,n);  %labels
ytest = zeros(m2,1);
%convert labels to -1 or 1
for i=1:m2
    if strcmp(yapp_sta{1}, ytest_RAW(i))
        ytest(i) = 1;
    else
        ytest(i) = -1;
    end;
end;

ypredapp = svmval_categorical(xtest,xsup,w,w0,kernel,kerneloption, symbolSta, dotProduct, my_lambda);
ypredapp = sign(ypredapp);    % Sign of ypredapp

correct_num = sum(ytest==ypredapp);
%disp('accuracy:');
correct_percent = correct_num / m2;  % Accuracy

%Compute Precision index
num1 = sum(ytest==1);  %number of positive class   %正类样本数
num2 = sum(ytest==-1); %number of negative  class  %负类样本数
TP1 = 0 ; % TP of positive class %正类的正确肯定
FP1 = 0 ; % FP of positive class %正类的错误肯定
FN1 = 0 ; % FN of positive class %正类的错误否定
TP2 = 0 ; % TP of negative class %负类的正确肯定
FP2 = 0 ; % FP of negative class %负类的错误肯定
FN2 = 0 ; % FN of negative class %负类的错误否定
for i=1:m2
    if(ytest(i)==1 && ypredapp(i)==1)  % The prediction is true and the actual is positive %预测为真，实际为真 
        TP1 = TP1 + 1;
    end;
    if(ytest(i)==-1 && ypredapp(i)==1) %The prediction is true and the actual is negative %预测为真，实际为假
        FP1 = FP1 + 1;
    end;
    if(ytest(i)==1 && ypredapp(i)==-1) %The prediction is false and the actual is positive %预测为假，实际为真
        FN1 = FN1 + 1;
    end;   
    
    
    if(ytest(i)==-1 && ypredapp(i)==-1)   
        TP2 = TP1 + 2;
    end;
    if(ytest(i)==1 && ypredapp(i)==-1) 
        FP2 = FP2 + 1;
    end;
    if(ytest(i)==-1 && ypredapp(i)==1) 
        FN2 = FN2 + 1;
    end;  
end;

Precision_1 = TP1/(TP1 + FP1);  % Precision of positive class %正类的查准率
Recall_1 = TP1/(TP1 + FN1); % Recall of positive class %正类的召回率
F1_Score_1 = 2*Precision_1* Recall_1/(Precision_1+Recall_1);

Precision_2 = TP2/(TP2 + FP2);   % Precision of negative class
Recall_2 = TP2/(TP2 + FN2);       %Recall of negative class
F1_Score_2 = 2*Precision_2* Recall_2/(Precision_2+Recall_2);

F1_Score = ( F1_Score_1*sum(ytest==1) + F1_Score_2*sum(ytest==-1) )/(sum(ytest==1)+sum(ytest==-1));
%fprintf('Weight F1_Measure: %6.3f\n', F1_Score);

result = F1_Score;
