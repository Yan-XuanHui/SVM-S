function [y]=svmval_categorical(x, xsup, w,b,kernel, kerneloption, dataSta, dotProduct,my_lambda)
% INPUT    
% x��test data set %Ҫ���Ե�����,�ж������
% xsup��support vector %֧������
% dataSta
% dotProduct

% ps: Kernel Matrix of x %�Ǻ˾���(����������֧��������xsup�ĺ˾���
ps = svmkernel_categorical(x, kernel, kerneloption, dataSta,dotProduct,my_lambda, xsup);
  
y = ps*w+b;
    

