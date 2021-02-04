function [M,V]=imlocalmv(I,r)
[m,n]=size(I);
x=1:m;
y=1:n;
[X,Y]=meshgrid(x,y);M=zeros(m,n);V=zeros(m,n);
for i=1:m
    for j=1:n
     M(i,j)= sum(sum(I(max(i-r,1):min(i+r,m),max(j-r,1):min(j+r,n))))/(min(i+r,m)-max(i-r,1)+1)/(min(j+r,n)-max(j-r,1)+1); 
     V(i,j)=sum(sum((I(max(i-r,1):min(i+r,m),max(j-r,1):min(j+r,n))-M(i,j)).^2))/(min(i+r,m)-max(i-r,1)+1)/(min(j+r,n)-max(j-r,1)+1);
    end
end
end