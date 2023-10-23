function [tHit,xHit,xHit_] = extractHits(simRes,varargin)
% extractHits - extracts time and states where guard intersection happened
%    default: all guard intersections, can be set to specific location
%    before the transition; note that we have to return cell-arrays for the
%    state vectors before and after jumping since subsequent locations need
%    not have the same number of states
%
% Syntax:
%    [tHit,xHit,xHit_] = extractHits(simRes)
%    [tHit,xHit,xHit_] = extractHits(simRes,locIDstart)
%    [tHit,xHit,xHit_] = extractHits(simRes,locIDstart,locIDend)
%
% Inputs:
%    simRes - simResult object
%    locIDstart - (optional) location vector where transition was triggered
%    locIDend - (optional) location vector to where state transitioned
%
% Outputs:
%    tHit - double vector of switching times (1xs)
%    xHit - cell-arrays of state vectors before switching (1xs)
%    xHit_ - cell-arrays of state vectors after switching (1xs)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simResult, hybridAutomaton/simulateRandom

% Authors:       Mark Wetzlinger
% Written:       21-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default values
[locIDstart,locIDend] = setDefaultValues({[],[]},varargin);

% check input arguments
inputArgsCheck({{simRes,'att','simResult'},...
                {locIDstart,'att','numeric'},...
                {locIDend,'att','numeric'}});

% empty object or all continuous-time simulation runs
if all(isempty(simRes)) ...
        || all(arrayfun(@(x) isscalar(x.loc) && x.loc == 0,simRes,'UniformOutput',true))
    tHit = []; xHit = {}; xHit_ = {}; return
end

tHit = zeros(0,0);
xHit = {};
xHit_ = {};

% loop over all trajectories
for r=1:length(simRes)

    % shortened code for getting all hits
    if isempty(locIDstart) && isempty(locIDend)
        tHit = [tHit; cellfun(@(x) x(1),simRes(r).t(2:end),'UniformOutput',true)];
        xHit = [xHit; cellfun(@(x) x(end,:)',simRes(r).x(1:end-1),'UniformOutput',false)];
        xHit_ = [xHit_; cellfun(@(x) x(1,:)',simRes(r).x(2:end),'UniformOutput',false)];
        continue
    end

    % longer version if locations have to be checked
    for j=1:length(simRes(r).t)-1
    
        % check if location before and after jump matches user input        
        startcond = isempty(locIDstart) || ...
            (all(size(simRes(r).loc(j,:)') == size(locIDstart)) && ...
            all(simRes(r).loc(j,:)' == locIDstart));
        endcond = isempty(locIDend) || ...
            (all(size(simRes(r).loc(j+1,:)') == size(locIDend)) && ...
            all(simRes(r).loc(j+1,:)' == locIDend));
    
        if startcond && endcond
            % time before and after jump are equal so doesn't matter
            tHit  = [tHit; simRes(r).t{j}(end)];
            % state vector before jump
            xHit = [xHit; simRes(r).x{j}(end,:)'];
            % state vector after jump
            xHit_ = [xHit_; simRes(r).x{j+1}(1,:)'];
        end    
    end

end

% ------------------------------ END OF CODE ------------------------------
