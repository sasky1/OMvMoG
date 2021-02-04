%无噪视频实验的视频数据读取
clc;clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%读取图片的数据，产生数据矩阵D
pathname = 'E:\学习\我的文件\RPCA\曹\Hall\';       
path = dir(pathname);
m=176;n=144;       %视频分辨率
c=150;             %视频帧数   
D = zeros(m*n,c);
for i=1:c
    str = strcat(pathname,path(i+3).name);
    A1 = imread(str);
    A = rgb2gray(A1);
    D(:,i) = A(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RPCA角度，调用算法IALM
[a,b]=size(D);
lambda=1/sqrt(a);
tic;[A_hat E_hat iter] = inexact_alm_rpca(D, lambda);toc;
D=uint8(D);A_hat=uint8(A_hat);E_hat=uint8(abs(E_hat));

%LRMF角度（随机初值）
[a,b]=size(D);InU=rand(a,1);D=double(D);tic;[OutU,OutV,t] = EfficientMF(D, InU,10, 0.00001);
A=OutU*OutV';E=abs(D-A);toc;



