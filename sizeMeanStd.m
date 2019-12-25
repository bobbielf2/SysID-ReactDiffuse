function U = sizeMeanStd(samps, threshold, dx, minmax)
% compute size distribution of samples

% samps = solution field to be sampled, ordered as "timestep-fast, species-slow" manner
%    e.g. say have n=2 species s1, s2, and k=3 snapshots i1, i2, i3
%    then order is samps = [s1(i1), s1(i2), s1(i3), s2(i1), s2(i2), s2(i3)]
% threshold = value of the level set that defines precipitates
% dx    = grid spacing

N = size(samps,1);  % total number of samples

% precompute min and max of samples
if nargin > 3 && size(minmax,2) == N
    cmin = minmax(1,:); cmax = minmax(2,:);
else
    cmin = min(samps,[],2); cmax = max(samps,[],2);
end

% compute size distribution
U = zeros(2,N);
for i = 1:N
    sample = samps(i,:);
%     if any(abs(sample)>1e2)
%         S = dx*numel(sample);
%         vals = [S;0]; % whole domain as one bump
%         U(:,i) = vals;
%     else
        % normalized sample to [0,1]
        sample_sc = (sample - cmin(i))/(cmax(i) - cmin(i));
        
        % find concentration sizes at a threshold
        A = (sample_sc >= threshold); % find continuous components
        A = diff(A); % starting & end points of components (starting labeled 1, ending labeled -1)
        startpt = find(A == 1); endpt = find(A == -1);
        if isempty(startpt) || isempty(endpt) || (numel(startpt)==1 && numel(endpt)==1 && endpt<startpt)
            S = dx*numel(sample);
            vals = [S;0]; % whole domain as one bump
            U(:,i) = vals;
        else
            try
                if endpt(1) < startpt(1), endpt(1) = []; end % remove incomplete component at left bdry
                if startpt(end) > endpt(end), startpt(end) = []; end % remove incomplete component  at right bdry
            catch
                keyboard
            end
            S = (endpt - startpt)*dx; % get sizes
            
            % size distribution
            vals = [mean(S);std(S)];
            U(:,i) = vals;
        end
%     end
end

end