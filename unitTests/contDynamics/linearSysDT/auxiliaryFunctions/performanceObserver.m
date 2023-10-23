function perf = performanceObserver(EstSet)
% performanceObserver - evaluates the performance of guaranteed
% state estimation of linear discrete-time systems according to [1].
%
% Syntax:
%    perf = performanceObserver(EstSet)
%
% Inputs:
%    EstSet - sets of estimated states
%
% Outputs:
%    perf - performance struct
%
% Reference:
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems, 
%        in preparation.

% Authors:       Matthias Althoff
% Written:       16-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Loop through all sets
% initialize values
perf.interval = [];
perf.max = [];
perf.min = [];
perf.vol = [];
perf.volBox = [];
perf.rad = [];
% no error propagation
for i = 1:length(EstSet.timePoint.set)
    try
        % compute interval hull
        perf.interval{end+1} = interval(EstSet.timePoint.set{i}); 
        % supremum of interval hull
        perf.max(:,end+1) = supremum(perf.interval{i}); 
        % infimum of interval hull
        perf.min(:,end+1) = infimum(perf.interval{i}); 
        % volume computation
        if isa(EstSet.timePoint.set{i},'zonotope')
            % approximate for zonotopes
            perf.vol(end+1) = volume(EstSet.timePoint.set{i},'alamo'); % approximate volume computation according to Alamo
        elseif isa(EstSet.timePoint.set{i},'ellipsoid')
            % exact for ellipsoids
            perf.vol(end+1) = volume(EstSet.timePoint.set{i});
        end
        % volume of interval hull
        perf.volBox(end+1) = volume(perf.interval{i}); 
        % radius of interval hull
        perf.rad(:,end+1) = rad(perf.interval{i});
    catch
        % interval probably too small
        disp('time point ignored for evaluation');
    end
end


% root mean square (rms) average of interval radii
for i = 1:length(perf.rad(:,1))
    perf.IRadius(i) = sqrt(mean(perf.rad(i,:).^2));
end

% ------------------------------ END OF CODE ------------------------------
