function [y]=svmval_categorical(x, xsup, w,b,kernel, kerneloption, dataSta, dotProduct,my_lambda)
% INPUT    
% x：test data set %要测试的数据,判断其类别
% xsup：support vector %支持向量
% dataSta
% dotProduct

% ps: Kernel Matrix of x %是核矩阵(测试数据与支持向量机xsup的核矩阵）
ps = svmkernel_categorical(x, kernel, kerneloption, dataSta,dotProduct,my_lambda, xsup);
  
y = ps*w+b;
    

