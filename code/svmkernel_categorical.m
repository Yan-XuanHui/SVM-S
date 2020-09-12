function [K,option]=svmkernel_categorical(x, kernel, kerneloption, dataSta, dotProduct, my_lambda, xsup)
% Computing Kernel between two objects and store in the matrix, %计算两两样本核运算结果，存在矩阵中

% Usage  K=svkernel_categorical(x,kernel,kerneloption);
%
% Returns the scalar product of the vectors x by using the
% mapping defined by the kernel function or x and xsup
% if the matrix xsup is defined
%
% Input
% 
% x		:input vectors
% kernel 	: kernel function
%		Type								Function					Option
%		Polynomial						'poly'					Degree (<x,xsup>+1)^d
%		Homogeneous polynomial		'polyhomog'				Degree <x,xsup>^d	
%		Gaussian							'gaussian'				Bandwidth
%		Heavy Tailed RBF				'htrbf'					[a,b]   %see Chappelle 1999	
%		Mexican 1D Wavelet 			'wavelet'
%		Frame kernel					'frame'					'sin','numerical'...	
%
%  kerneloption	: scalar or vector containing the option for the kernel
% 'gaussian' : scalar gamma is identical for all coordinates
%              otherwise is a vector of length equal to the number of 
%              coordinate
%
% 'poly' : kerneloption is a scalar given the degree of the polynomial
%          or is a vector which first element is the degree of the polynomial
%           and other elements gives the bandwidth of each dimension.
%          thus the vector is of size n+1 where n is the dimension of the problem.
%
%	see also svmreg,svmclass,svmval, kernelwavelet,kernelframe
%  dataSta: the matrix of data frequency statistic %符号频率统计矩阵
%  dotProduct: dot product matrix    %点积矩阵

k_lambda = my_lambda;  %band width 
if nargin<3
    kerneloption=1;
end;
if nargin<2
    kernel='gaussian';
end;
if nargin<7
    xsup=x;
end;
[n1 n2]=size(x);
[n n3]=size(xsup);
ps  =  zeros(n1,n);			% product scalaire
switch lower(kernel)        
case 'poly'
    
     [nk,nk2]=size(kerneloption);   
     if nk>nk2
         kerneloption=kerneloption';
         nk2=nk;
     end;
     if nk2==1
         degree=kerneloption;
         var=ones(1,n2);
         
     elseif nk2 == 2  
         degree=kerneloption(1);
         var=ones(1,n2)*kerneloption(2); % the second parameter of kerneloption %kerneloption的第2个参数命名为var
         
     elseif nk2== n2+1
         degree=kerneloption(1);
         var=kerneloption(2:n2+1);
         
     elseif nk2 ==n2+2
         degree=kerneloption(1);
         var=kerneloption(2:n2+1);
     end;
 
     if nk2==1
         aux=1;
     else
         aux=repmat(var,n,1);
     end;
   
%     ps= x *(xsup.*aux.^2)';

    % if there is parameter xsup, computing dotProduct  %如果有参数xsup,重新计算dotProcuct
    if (nargin>6)
       % computing dot product
       [dotProduct] = dotProductMatrix(x, dataSta, my_lambda, xsup); 
    end;
    ps = dotProduct;
    if degree > 1
        K = (ps+1).^degree;
    else
        K = ps;
    end;
case 'polyhomog'
    
    [nk,nk2]=size(kerneloption);   
    if nk>nk2
        kerneloption=kerneloption';
        nk2=nk;
    end;
    if nk2==1
        degree=kerneloption;
        var=ones(1,n2);
    else
        if nk2 ~=n2+1
            degree=kerneloption(1);
            var=ones(1,n2)*kerneloption(2);
        else
            degree=kerneloption(1);
            var=kerneloption(2:nk2);
        end;
    end;
    
    aux=repmat(var,n,1);
    ps= x *(xsup.*aux.^2)';
    K =(ps).^degree;
        
case 'gaussian'
%     [nk,nk2]=size(kerneloption);
%     if nk ~=nk2
%         if nk>nk2
%             kerneloption=kerneloption';
%         end;
%     else
%         kerneloption=ones(1,n2)*kerneloption;
%     end;
%     
%     if length(kerneloption)~=n2 & length(kerneloption)~=n2+1 
%         error('Number of kerneloption is not compatible with data...');
%     end;
    
    % computing Gaussian Kernel %计算高斯核函数的值
    [N, dim] = size(x);
    [M, dim1] = size(xsup);
    ps = zeros(N, M);
	
	distance_type = 1;    %distance_type = 1: sqrt(2*(1-lambda)^2 )   distance_type = 2: frequency difference forumal %用频率差异公式(用在稀有类挖掘中)
    %计算ps的每个元素值
    for i=1:N      %遍历x
        for j=1:M  %遍历xsup
            ds = zeros(1,dim);     
            for d=1:dim   
                xd = x(i,d);        %第d维对应的符号
                yd = xsup(j,d);
                if strcmp(xd, yd)   %相等直接距离为0
                    ds(d) = 0;
                else       
					if(distance_type==1)
						ds(d) =  sqrt(2*(1-k_lambda(d))^2 );
					else								
						% computing the frequcncy of xd  %计算xd符号的频率
						[bool, x_idx] = ismember(xd, dataSta{d}(:,1));
						fx = dataSta{d}(:,3);      % the 3th colmun of dataSta is frequency*100  %dataSta的第3列是频率*100
						fxd = fx{x_idx} / 100.0;   % frequency of xd  %xd的频率
						% computing the frequcncy of yd
						[bool, y_idx] = ismember(yd, dataSta{d}(:,1));
						fy = dataSta{d}(:,3);
						fyd = fy{y_idx} / 100.0;   % frequency of yd 
						ds(d) = (abs(fxd - fyd) +  (1.0/2/N))*(1-k_lambda(d)) ; % using frequcncy difference formula %利用频率差距离公式
					end; %if(distance_type==1)
                end;%if strcpm(xd,yd)
                
            end; %for d
            ps(i, j) = sum( ds.^2)/(kerneloption^2);
        end; %for j
    end; %for i

    K = exp(-ps/2);

case 'htrbf'    % heavy tailed RBF  %see Chappelle Paper%
    b=kerneloption(2);
    a=kerneloption(1);
    for i=1:n
        ps(:,i) = sum( abs((x.^a - ones(n1,1)*xsup(i,:).^a)).^b   ,2);
    end;   
    K = exp(-ps);
    
case 'gaussianslow'    %
    %b=kerneloption(2);
    %a=kerneloption(1);
    for i=1:n
        ps(:,i) = sum( abs((x - ones(n1,1)*xsup(i,:))).^2 ,2)./kerneloption.^2/2;
    end;
    K = exp(-ps);
case 'multiquadric'
    metric = diag(1./kerneloption);
    ps = x*metric*xsup'; 
    [nps,pps]=size(ps);
    normx = sum(x.^2*metric,2);
    normxsup = sum(xsup.^2*metric,2);
    ps = -2*ps + repmat(normx,1,pps) + repmat(normxsup',nps,1) ; 
    K=sqrt(ps + 0.1);
case 'wavelet'
    K=kernelwavelet(x,kerneloption,xsup);     
case 'frame'
    K=kernelframe(x,kerneloption,xsup,framematrix,vector,dual);
case 'wavelet2d'
    K=wav2dkernelint(x,xsup,kerneloption);
case 'radialwavelet2d'
    K=radialwavkernel(x,xsup);    
case 'tensorwavkernel'
    [K,option]=tensorwavkernel(x,xsup,kerneloption);  

case 'numerical'
    K=kerneloption.matrix;
case 'polymetric'
    K=x*kerneloption.metric*xsup';
    
case 'jcb'
    K=x*xsup';
    
end;