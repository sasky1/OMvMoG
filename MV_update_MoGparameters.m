function [model,label,R,h,v,vv]=MV_update_MoGparameters(x,model,lamda,w_x,vv,h)
V=size(x,2);
  for i=1:V  
U{i}=model{i}.U;
N{i}=model{i}.N;
nmodel{i}.weight=model{i}.weight;
nmodel{i}.mu=model{i}.mu;
nmodel{i}.Sigma=model{i}.Sigma;
  end
r=size(U{1},2);
if nargin <4
for i=1:V
w_x{i}=ones(size(x{i}));
end
end
for i=1:V
ind{i}= (w_x{i}==0);
m{i}=size(x{i},1);
h{i}=ones(m{i},1);h{i}(ind{i})=0;
end

iter=10;
%     vv=zeros(r,1);
for tt=1:iter
  for i=1:V  
v{i}=(U{i}'.*repmat(h{i}',r,1)*U{i}+lamda*eye(r))^-1*(U{i}'*(h{i}.*x{i})+lamda*vv);
%E_step
[R{i}] = expectation((x{i}-U{i}*v{i})', nmodel{i});
[~,label{i}]=max(R{i}');
label{i}=label{i}';
%M_step
%M_step for w sigma N
[nmodel{i}] = maximizationModel((x{i}-U{i}*v{i})',R{i},model{i},N{i});
h{i}=sum(R{i}./repmat(2*nmodel{i}.Sigma,m{i},1),2);h{i}(ind{i})=0;
% v=(U'.*repmat(h',r,1)*U+0.0001*eye(r))^-1*U'*(h.*x);
  end
  vv=zeros(r,1);
  for i=1:V
  vv=vv+v{i};
  end
  vv=vv/V;
end
  for i=1:V  
model{i}.weight=nmodel{i}.weight;
model{i}.mu=nmodel{i}.mu;
model{i}.Sigma=nmodel{i}.Sigma;
  end
end

function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;
end

function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
logRho = zeros(n,k);


for i = 1:k
    logRho(:,i) = loggausspdf(X,mu(i),Sigma(i));
end
logRho = bsxfun(@plus,logRho,log(w));
T = logsumexp(logRho,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);
end

function [nmodel] = maximizationModel(X,R,model,N)
forgetnumber=0;
ww=model.weight;
 mmu=model.mu;
ssigma=model.Sigma;

[d,n] = size(X);
k = size(R,2);

nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
% mu = bsxfun(@times, X*R, 1./nk);
w = nk/size(R,1);
Sigma = zeros(1,k);
sqrtR = sqrt(R);
new_Nk=(1-forgetnumber)*ww*N+(1+forgetnumber)*sum(R);
new_N=(1-forgetnumber)*N+(1+forgetnumber)*size(X,2);
new_w=new_Nk/new_N;
for i = 1:k
    Xo = bsxfun(@minus,X,mu(i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(i) = ((1+forgetnumber).^2*Xo*Xo'+(1-forgetnumber).^2*ssigma(i)*ww(i)*N)/new_Nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end


nmodel.mu = mu;
nmodel.Sigma =Sigma;
nmodel.weight = new_w;
end



