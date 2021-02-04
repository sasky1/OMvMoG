function [label,mu,sigma,llk] = maxguass(data,k)
%UNTITLED23 Summary of this function goes here
%   Detailed explanation goes here


[m,n]=size(data);

%   [~,mu]= kmeans(data',k);
% for i=1:k
%     dd=randperm(n);
%   mu(i,:)= sum(data(:,dd(1:fix(n/k)))')/fix(n/k);
% 
% end
  mu=zeros(1,k);
sigma=0.0005*rand(1,k)+0.001;
% for i=1:k
% sigma(i)=1;
% end

iter=100;eps=0.000001;llk=0;
for t=1:iter

p=[];
for i=1:k
  p=[p, log(mvnpdf(data',mu(i),(sigma(i)+sigma(i)')/2))];
end



[~,label]=max(p');




for i=1:k
%    mu(i)= sum(data(:,find(label==i))')/length(find(label==i));
  d=data(:,find(label==i))'-repmat(mu(i),length(find(label==i)),1);
  sigma(i)=d'*d/length(find(label==i));
  sigma(i)=(sigma(i)+sigma(i)')/2;
end

llk0=llk;
llk=sum(max(p'));
if abs(llk-llk0)/abs(llk)<eps;
   break; 
end


end
llk=llk/n;
