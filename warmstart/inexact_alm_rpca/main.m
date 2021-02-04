%������Ƶʵ�����Ƶ���ݶ�ȡ
clc;clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ȡͼƬ�����ݣ��������ݾ���D
pathname = 'E:\ѧϰ\�ҵ��ļ�\RPCA\��\Hall\';       
path = dir(pathname);
m=176;n=144;       %��Ƶ�ֱ���
c=150;             %��Ƶ֡��   
D = zeros(m*n,c);
for i=1:c
    str = strcat(pathname,path(i+3).name);
    A1 = imread(str);
    A = rgb2gray(A1);
    D(:,i) = A(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RPCA�Ƕȣ������㷨IALM
[a,b]=size(D);
lambda=1/sqrt(a);
tic;[A_hat E_hat iter] = inexact_alm_rpca(D, lambda);toc;
D=uint8(D);A_hat=uint8(A_hat);E_hat=uint8(abs(E_hat));

%LRMF�Ƕȣ������ֵ��
[a,b]=size(D);InU=rand(a,1);D=double(D);tic;[OutU,OutV,t] = EfficientMF(D, InU,10, 0.00001);
A=OutU*OutV';E=abs(D-A);toc;



