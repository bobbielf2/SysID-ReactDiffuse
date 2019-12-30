function val = my_lnLikelihood(theta, data)

% Note the data variable is 4 x number of snapshots * number of species.
ydata = data.ydata;
full_QoI_sigmas = data.full_QoI_sigmas;

pred = data.G(theta);

% Uses additive independent Gaussian noise with fixed sigmas.
val = 0.0;
for i = 1:size(ydata, 1)
  % Loops through the individual QoIs.
  
  for j = 1:size(ydata, 2)
    % Loops through the snapshots.
  
    val = val - 0.5 * log(2 * pi) - log( full_QoI_sigmas(i,j) ) - 0.5 * ( ...
        (ydata(i,j)-pred(i,j)) / full_QoI_sigmas(i,j) ).^2; 
    
  end
  
end
