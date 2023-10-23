function out = masks(phi,dom,val)
% masks - calculate signals describing the time intervals in which the
% value of each atomic proposition affects the valuation of the formula.
%
% The signals have as duration the horizon of the formula.
%
% Syntax:
%    m = masks(phi);
%
% Inputs:
%    phi - STL formula
%    dom - time-domain interval we want to monitor
%    val - value the atomic proposition is supposed to have
%
% Outputs:
%    out - map of atomic propositions to signals
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       23-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

arguments
    phi stl
    dom interval = interval(0,0)
    val string = "any"
end

map = aux_collectIntervals(phi, dom, true);
k = keys(map);
out = containers.Map;

dur = supremum(dom + maximumTime(phi));

% combine all intervals to a single signal
for i = 1:length(k)
    c = map(k{i});
    ints = c{1};
    dirs = c{2};

    out(k{i}) = signal(dur, false);

    for j = 1:length(ints)
        if val == "any" || val == string(dirs(j))
            out(k{i}) = out(k{i}) | signal.indicator(dur, ints(j), true);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function out = aux_collectIntervals(phi, dom, dir)
    % Collect all intervals for every atomic proposition where the value of
    % that proposition matters for the evaluation of the formula.

    if strcmp(phi.type, 'true')
        out = containers.Map;
    elseif strcmp(phi.type, 'false')
        out = containers.Map;
    elseif strcmp(phi.type, 'variable')
        out = containers.Map(formattedDisplayText(phi), {dom dir});
    elseif strcmp(phi.type, '&')
        lhs = aux_collectIntervals(phi.lhs, dom, dir);
        rhs = aux_collectIntervals(phi.rhs, dom, dir);

        out = aux_combineMaps(lhs, rhs);
    elseif strcmp(phi.type, '|')
        lhs = aux_collectIntervals(phi.lhs, dom, dir);
        rhs = aux_collectIntervals(phi.rhs, dom, dir);

        out = aux_combineMaps(lhs, rhs);
    elseif strcmp(phi.type, '~')
        out = aux_collectIntervals(phi.lhs, dom, ~dir);
    elseif strcmp(phi.type, 'until')
        m1 = aux_collectIntervals(phi.lhs, dom + interval(0, phi.to), dir);
        m2 = aux_collectIntervals(phi.rhs, dom + interval(phi.from, phi.to), dir);

        out = aux_combineMaps(m1, m2);
    elseif strcmp(phi.type, 'release')
        m1 = aux_collectIntervals(phi.lhs, dom + interval(0, phi.to), dir);
        m2 = aux_collectIntervals(phi.rhs, dom + interval(phi.from, phi.to), dir);

        out = aux_combineMaps(m1, m2);
    elseif strcmp(phi.type, 'finally')
        out = aux_collectIntervals(phi.lhs, dom + interval(phi.from, phi.to), dir);
    elseif strcmp(phi.type, 'globally')
        out = aux_collectIntervals(phi.lhs, dom + interval(phi.from, phi.to), dir);
    end
end

function m1 = aux_combineMaps(m1, m2)
    % Combine two maps by combining all values

    k = keys(m2);

    for i = 1:length(k)
        if isKey(m1, k{i})
            m1(k{i}) = cellfun(@(a,b) [a b], m1(k{i}), m2(k{i}), 'Uniform', 0);
        else
            m1(k{i}) = m2(k{i});
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
