classdef simResultAnalyzer < handle
% simResultAnalyzer - construct boolean signals for the given atomic
% propositions over simulation results
%
% Syntax:
%    analyzer = simResultAnalyzer(aps, 10)
%
% Inputs:
%    aps - array of sets representing the atomic propositions
%    dur - duration of the signal
%    index - index of the simulation to analyze
%
% Outputs:
%    obj - simResultAnalyzer object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Benedikt Seidl
% Written:       04-January-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % map of propositions to their geometric representation
    aps containers.Map

    % the simulation result to analyze
    simRes
end

methods

    function obj = simResultAnalyzer(aps, simRes)
        obj.aps = aps;
        obj.simRes = simRes;
    end

    function sigs = analyze(obj, n)
        % analyze all time intervals of simulation result with index n
        k = keys(obj.aps);
        sigs = containers.Map;
        dur = obj.simRes.t{n}(end);

        % prepare signals for all atomic propositions
        for i = 1:length(k)
            sigs(k{i}) = signal(dur, false);
        end

        for i = 1:length(obj.simRes.t{n})-1
            % obtain value of simulation
            %
            % Taking the first point only is an approximation.
            % However, since simulation results have no claim of
            % correctness this inaccuracy is acceptable here.
            point = obj.simRes.x{n}(i,:);

            % set location property
            if ~isempty(obj.simRes.loc)
                loc = obj.simRes.loc(n);
            else
                loc = 1;
            end

            % calculate time interval
            from = obj.simRes.t{n}(i);
            to = obj.simRes.t{n}(i+1);

            % check all atomic propositions
            for j = 1:length(k)
                if obj.aps(k{j}).evaluatePoint(point.', loc)
                    sig = signal.indicator(dur, interval(from, to), true);

                    sigs(k{j}) = sigs(k{j}) | sig;
                end
            end
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
