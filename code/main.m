
% Description: The main program of SVM-S alogrithm test
% Author: Xuanhui Yan
% Data: 2019-04-01

clear;
run_count = 5; %number of test
avg_result = linspace(1, run_count, run_count);
class = 2 ;   % The number of class in the data set 
filePath = '..\data\BreastCancer.xlsx';    % file path of data set 
C = 10 ;      %parameter C of SVM
i = 1;
while i <= run_count
    
    fprintf('%dst test:\n',i);
    if(class<=2)
        %result = SVM_Categorical_Fun(filePath,'gaussian', 0.5, C);    % using gaussian Kernel
        result = SVM_Categorical_Fun(filePath, 'poly', 3, C);          % using poly Kernel
    else
        %result = SVM_Categorical_MultiClass_Fun(filePath, 'gaussian', 0.6, C);
        result = SVM_Categorical_MultiClass_Fun(filePath, 'poly', 3, C);
    end;
    avg_result(i) = result;
    fprintf('WF1 = %f\n', result);
    i = i+1;
end
avg_result
fprintf('Avg. WF1_Score=%f\n', mean(avg_result));