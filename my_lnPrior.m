function val = my_lnPrior(theta, params)

% Gaussian prior here.
mu0 = params.parmu0;                       % prior means
sig0  = params.parsig0;                      % prior stds

m = length(theta);  % Parameter dimension.
%val = - m/2 * log(2 * pi) - sum( log(sig0) ) - 0.5 * sum( (theta-mu0 ./ sig0).^2 );
val = 1;