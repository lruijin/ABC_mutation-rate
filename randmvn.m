function D = randmvn(mu,sigma,N,eps)
%RANDMVN Random matrices from the multivariate normal distribution.
%   D = RANDMVN(MU,SIGMA,N) returns a matrix of random samples chosen   
%   from the multivariate normal distribution with mean vector, MU, and 
%   covariance matrix, SIGMA. N is the number of samples (columns in D).
%
%   SIGMA is a symmetric positive definite matrix with size equal to the 
%   length of MU
%%%%%%%%%%%%%%%%%%%%%%% MU must be a column vector %%%%%%%%%%%%%%%%%%%
%[m1 n1] = size(mu);
%c = max([m1 n1]);
%if m1 .* n1 ~= c
%   error('MU must be a vector.');
%end
%if n1 == c
%  mu = mu';
%end
%%%%%%%%%%%%%%%%%%%%%%% SIGMA must be a cov. matrix %%%%%%%%%%%%%%%%%%
%[m n] = size(sigma);
%if m ~= n
%   error('SIGMA must be square');
%end

%if m ~= c
%   error('The length of MU must equal the number of rows in SIGMA.');
%end

%[T p] = chol(sigma);
[T,p] = jitterChol(sigma,eps);
if p ~= 0
  sigma_new = diag(diag(sigma));
  [T,~] = jitterChol(sigma_new,eps);
  %error('SIGMA must be a positive definite matrix.');
end
%%%%%%%%%%%%%%%%%%%%%%% calculating D %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = mu(:,ones(N,1));  % repeating mu in each column

%D = randn(N,c) * T + mu;
lam = size(T,1);
if lam ~= size(mu, 1)
    error(['size of T is', num2str(lam)]);
end
D = (randn(N,size(mu,1)) * T)' + mu;
