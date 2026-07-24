function C = mwdct(ntpts,K)

% function C = dct(ntpts,K)
%
% Creates Discrete Cosine Transform basis set
%
% ntpts - number of timepoints
% K - number of cosines

n = (0:(ntpts-1))';

C = zeros(size(n,1),K);

C(:,1)=ones(size(n,1),1)/sqrt(ntpts);
for k=2:K
  C(:,k) = sqrt(2/ntpts)*cos(pi*(2*n+1)*(k-1)/(2*ntpts));
end

