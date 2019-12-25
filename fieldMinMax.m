function U = fieldMinMax(samps)
% compute matrix composition (mean of local minima) of samples

% samps = solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%    e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%    then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
% qoi   = struct containing parameters needed to compute the QoI. (See qoiInit.m)

N = size(samps,1);  % total number of samples

% compute matrix composition (mean of local minima)
U = zeros(2,N); 
for i = 1:N
    sample = samps(i,:);
    % find local max & min
    variation = sign(diff(sample));
    extrema = sign(diff(variation));
    imax = find(extrema == -1) + 1;
    imin = find(extrema == 1) + 1;
    
    if ~isempty(imin)
        U(1,i) = mean(sample(imin)); % mean of loc mins
    else
        U(1,i) = min(sample,[],'all');
    end
    if ~isempty(imax)
        U(2,i) = mean(sample(imax)); % mean of loc mins
    else
        U(2,i) = max(sample,[],'all');
    end
    % handle anomaly
    if isnan(U(1,i)), U(1,i) = -1e5; end
    if isnan(U(2,i)), U(2,i) = 1e5; end
end

end