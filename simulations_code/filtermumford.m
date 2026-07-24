function [F]=filtermumford(cut,TR,Nscans)


 %Construct the HP filter.  This should exactly match FSL
 % cut is the cutoff IN SECONDS that would be specified in FSL
 
 %cut=200;
 
 %TR is the TR in seconds
 %TR=3;
 
 % Nscans is the number of TRs (# of rows in design matrix).
 %Nscans=200;
 
 cut=cut/TR;
 sigN2=(cut/sqrt(2))^2; %variance approximation (???)
  K    = toeplitz(1/sqrt(2*pi*sigN2)*exp(-[0:(Nscans - 1)].^2/(2*sigN2))); %create a #TR x #TRs toeplitz matrix using the Gaussian distribution pdf
  K    = spdiags(1./sum(K')',0,Nscans,Nscans)*K;   %create a #TR x #TR sparse matrix (with 1/sum(K')':this create an inverse Gaussian pdf) and multiply by the toeplitz matrix

  
  % Smoothing matrix, s.t. H*y is smooth line
  H = zeros(Nscans,Nscans); 
  X = [ones(Nscans,1) (1:Nscans)'];
  
  for k = 1:Nscans
  W = diag(K(k,:));
  Hat = X*pinv(W*X)*W;
  H(k,:) = Hat(k,:);
  end
  
 %F is the filtering matrix that you premultiply the data and design by
 F=eye(Nscans)-H;