classdef trajectoryAnalyzer < handle
% trajectoryAnalyzer - construct boolean signals for the given atomic
% propositions over simulation results
%
% Syntax:
%    analyzer = trajectoryAnalyzer(aps, 10)
%
% Inputs:
%    aps - array of sets representing the atomic propositions
%    dur - duration of the signal
%    s - index of the trajectory realization to analyze
%
% Outputs:
%    obj - trajectoryAnalyzer object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Benedikt Seidl, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % map of propositions to their geometric representation
    aps containers.Map

    % the simulation result to analyze
    traj
end

methods

    function obj = trajectoryAnalyzer(aps, traj)
        obj.aps = aps;
        obj.traj = traj;
    end

    function sigs = analyze(obj, s, j)
        if nargin > 2 && length(obj) >= j
            obj = obj(j);
        end

        % analyze all time intervals of trajectory realization s
        k = keys(obj.aps);
        sigs = containers.Map;
        dur = obj.traj.t(:,end);

        % prepare signals for all atomic propositions
        for m = 1:length(k)
            sigs(k{m}) = finiteSignal(dur, false);
        end

        for l = 1:size(obj.traj.t,2)-1
            % obtain value of simulation
            %
            % Taking the first point only is an approximation.
            % However, since simulation results have no claim of
            % correctness this inaccuracy is acceptable here.
            point = obj.traj.x(:,l,s);

            % set location property
            if ~isempty(obj.traj.loc)
                loc = obj.traj.loc(:,l);
            else
                loc = 1;
            end

            % calculate time interval
            from = obj.traj.t(:,l);
            to = obj.traj.t(:,l+1);

            % check all atomic propositions
            for m = 1:length(k)
                if obj.aps(k{m}).evaluatePoint(point, loc)
                    sig = finiteSignal.indicator(dur, interval(from, to), true);
                    sigs(k{m}) = sigs(k{m}) | sig;
                end
            end
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
