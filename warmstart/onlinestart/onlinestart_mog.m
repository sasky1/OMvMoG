function [label,model,W,OutU,OutV,llh,A,B] = onlinestart(InW,InX,r,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deyu Meng, Fernando De la Torre. Robust matrix factorization %
% with unknown noise, ICCV, 2013                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform EM algorithm for fitting the MoG model.
% Step 1: Max parameters;
% Step 2: Expectation;
% Step 3: Max weighted L2 MF;
% Step 4: Expectation.
%USAGE:[label,model,W,OutU,OutV,llh] = MLGMDN(InW,InX,r,param)
%Input:
%   InX: d x n input data matrix
%   InW: d x n indicator matrix
%   r:   the rank
%   param.method: 1 Damped Newton method; 2 CoD method; 3 ALM method;
%                 4 Weighted PCA method.
%   param.maxiter: maximal iteration number
%   param.OriX: ground truth matrix
%   param.InU,InV: Initialized factorized matrices
%   param.k: the number of GMM
%   param.display: display the iterative process
%   param.tol: the thresholding for stop
%Output:
%   label:the labels of the noises
%   model:model.mu, the means of the different Gaussians
%         model.sigma,the variance of the different Gaussians
%         model.weight,the mixing coefficients
%   W: d x n weighted matrix
%   OutU: the fianl factorized matrix U
%   OutV: the fianl factorized matrix V
%   llh:  the log likelihood
%
% Author: Deyu Meng(dymeng@mail.xjtu.edu.cn) and Hongwei Yong
% Date: 10-31-2014


%% initialization
[m n] = size(InX);
if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end

if (~isfield(param,'OriX'))
    OriX = InX;
else
    OriX = param.OriX;
    clear param.OriX;
end

IND = find(InW(:) ~= 0);
if (~isfield(param,'InU'))
    s = median(abs(InX(IND)));
    s = sqrt(s/r);
    if min(InX(IND)) >= 0
        InU = rand(m,r)*s;
    else
        InU = rand(m,r)*s*2-s;
    end
else
    InU = param.InU;
end

if (~isfield(param,'InV'))
    if min(InX(IND)) >= 0
        InV = rand(n,r)*s;
    else
        InV = rand(n,r)*s*2-s;
    end
else
    InV = param.InV;
end

if (~isfield(param,'k'))
    k = 3;
else
    k = param.k;
end

if (~isfield(param,'display'))
    display = 0;
else
    display = param.display;
end

if (~isfield(param,'NumIter'))
    NumIter = 40;
else
    NumIter = param.NumIter;
end

if (~isfield(param,'tol'))
    tol = 1.0e-7;
else
    tol = param.tol;
end


%%Initialize GMM parameters
IND = find(InW(:) ~= 0);
tempX=InX(IND);
R = initialization(tempX',k);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));
model.mu = zeros(1,k);
model.Sigma = rand(1,k);
nk = sum(R,1);
model.weight = nk/size(R,1);
% llh = -inf(1,maxiter);
converged = false;
TempU = InU;
TempV = InV;
TempX = TempU*TempV';
Error = InX(:)-TempX(:);
Error = Error(IND);

t = 1;
%%%%%%%%%%%%%%%%Initialized E Step %%%%%%%%%%%%%%%%%%%
[R, llh(t)] = expectation(Error',model);
%%%%%%%%%%%%%%%%Initialized E Step %%%%%%%%%%%%%%%%%%%

while ~converged && t < maxiter
    t = t+1;
    
    %%%%%%%%%%%%%%%% M Step 1 %%%%%%%%%%%%%%%%%%%
    [model] = maximizationModel(Error',R);
    %%%%%%%%%%%%%%%% M Step 1 %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    [R, llh(t)] = expectation(Error',model);
    L1 = llh(t);
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% M Step 2 %%%%%%%%%%%%%%%%%%%
    [W TempU TempV ] = maximizationW(model,InW,InX,TempU,TempV,r,R,NumIter,param);
    TempX = TempU*TempV';
    Error = InX(:)-TempX(:);
    Error = Error(IND);
    %%%%%%%%%%%%%%%% M Step 2 %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    [R, llh(t)] = expectation(Error',model);
    L2 = llh(t);
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    
    if display
        disp(['The likelihood in this step is ',num2str(L1),' and ',num2str(L2),';']);
        disp(['There are ',num2str(k),' Gaussian noises mixed in data']);
        disp(['with variances ',num2str(model.Sigma)]);
        disp(['with weights ',num2str(model.weight),'.']);
        if (~isfield(param,'orimissing'))
            disp(['Relative reconstruction error ', num2str(sum(sum(((OriX - TempU*TempV')).^2))/sum(sum((OriX).^2)))]);
        else
            disp(['Relative reconstruction error ', num2str(sum(sum((InW.*(OriX - TempU*TempV')).^2))/sum(sum((InW.*OriX).^2)))]);
        end
        disp(['L2 RMSE is ', num2str(sqrt(mean(Error.^2)))]);
        disp(['L1 RMSE is ', num2str(mean(abs(Error)))]);
    end
    
    if (~isfield(param,'k'))
        [~,label(:)] = max(R,[],2);
        u = unique(label);   % non-empty components
        if size(R,2) ~= size(u,2)
            R = R(:,u);   % remove empty components
        else
            converged = llh(t)-llh(t-1) < tol*abs(llh(t));
        end
        k = length(u);
    else
        converged = llh(t)-llh(t-1) < tol*abs(llh(t));
    end
    
end
OutU = TempU;
OutV = TempV;

A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(OutV'*diag(W(i,:))*OutV)^-1;
    B(:,i)=OutV'*diag(W(i,:))*(InX(i,:))';
end


if ~display
    disp(['The likelihood in this step is ',num2str(L1),' and ',num2str(L2),';']);
    disp(['There are ',num2str(k),' Gaussian noises mixed in data']);
    disp(['with variances ',num2str(model.Sigma)]);
    disp(['with weights ',num2str(model.weight),'.']);
    if (~isfield(param,'orimissing'))
        disp(['Relative reconstruction error ', num2str(sum(sum(((OriX - TempU*TempV')).^2))/sum(sum((OriX).^2)))]);
    else
        disp(['Relative reconstruction error ', num2str(sum(sum((InW.*(OriX - TempU*TempV')).^2))/sum(sum((InW.*OriX).^2)))]);
    end
    disp(['L2 RMSE is ', num2str(sqrt(mean(Error.^2)))]);
    disp(['L1 RMSE is ', num2str(mean(abs(Error)))]);
end
if converged
    fprintf('Converged in %d steps.\n',t-1);
else
    fprintf('Not converged in %d steps.\n',maxiter);
end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
    while k ~= length(u)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
        [u,~,label] = unique(label);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d  %initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
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

function [model] = maximizationModel(X,R)

[d,n] = size(X);
k = size(R,2);

nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
% mu = bsxfun(@times, X*R, 1./nk);
w = nk/size(R,1);
Sigma = zeros(1,k);
sqrtR = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(i) = Xo*Xo'/nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;


function [W TempU TempV] = maximizationW(model,InW,InX,TempU,TempV,r,R,NumIter,param)
IND = find(InW(:)~=0);
k = size(R,2);
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;
r = size(TempU,2);
W = zeros(size(InX));
C = zeros(size(InX));
for j = 1:k
    W(IND) = W(IND) + R(:,j)/(2*Sigma(j));
    C(IND) = C(IND) + R(:,j)*mu(j)/(2*Sigma(j));
end
C(IND) = C(IND)./W(IND);
[TempU,TempV] = EfficientMCL2(InX-C, sqrt(W), TempU,TempV, NumIter, 0.00000001);
%     [TempU,TempV] = ALMwmf(InX-C, sqrt(W), TempU,TempV, NumIter, 0.00000001);






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
