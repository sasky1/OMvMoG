function model=MV_update_subspace(x,h,v,model)
V=size(x,2);
for i=1:V
U=model{i}.U;
A=model{i}.A;
B=model{i}.B;
ro=model{i}.lamda;
temp=sum(A.*repmat(v{i}',size(A,1),1,size(A,3)),2)/ro;
temp1=sum(A.*repmat(v{i},1,size(A,2),size(A,3)),1)/ro;
new_A=A/ro-repmat(reshape(h{i},1,1,size(h{i},1)),size(A,1),size(A,2),1).*repmat(temp,1,size(A,2),1).*repmat(temp1,size(A,1),1,1) ...
 ./(1+repmat(reshape(h{i},1,1,size(h{i},1)),size(A,1),size(A,2),1).*repmat(sum(temp.*repmat(v{i},1,size(temp,2),size(temp,3)),1),size(A,1),size(A,2),1));   
 new_B=B*ro+repmat((h{i}.*x{i})',size(B,1),1).*repmat(v{i},1,size(B,2));
model{i}.U=(reshape(sum(A.*repmat(reshape(B,1,size(A,2),size(A,3)),size(A,1),1,1),2),size(U,2),size(U,1),1))'; 
%output model.w model.sigma U v A,B,N
model{i}.A=new_A;
model{i}.B=new_B;
end

