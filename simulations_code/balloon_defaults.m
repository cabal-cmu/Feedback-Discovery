function [balloon_priors]=balloon_defaults()

% kappa    gamma   
balloon_priors.balloon_cbf_mean=[0.65, 0.41]';
% tau    alpha    E0  
balloon_priors.balloon_mean=[0.98, 0.32, 0.34]';
% logepsilon
balloon_priors.balloon2_mean=[0.00]';
    
balloon_priors.balloon_cbf_cov=[0.0148  -0.0001; -0.0001 0.0018];
balloon_priors.balloon_cov=[0.0514  0.0021  0.0018;  0.0021  0.0002  0.0001; 0.0018  0.0001  0.0001];
balloon_priors.balloon2_cov=[0.0312];
