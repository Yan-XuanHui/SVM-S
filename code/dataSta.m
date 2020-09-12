function [dataSta] = dataSta(xapp)

[r, dim] = size(xapp);		% r£ºrows£¬dim£ºcolmuns

%%
%% function of data frequency statistic for ecah attribute£¨using tabulate.m in statistic toolbox£©
for i=1:dim 
    tt = xapp(:,i);
    B{1,i} = tabulate(tt);
end;

dataSta = B;