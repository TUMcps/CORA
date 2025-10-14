function [tHit,xHit,xHit_,idzHit] = extractHits(traj,varargin)
% extractHits - extracts time and states where guard intersection happened
%    (i.e., location changes)
%    default: all guard intersections, can be set to specific location
%    before the transition
%
% Syntax:
%    [tHit,xHit,xHit_] = extractHits(traj)
%    [tHit,xHit,xHit_] = extractHits(traj,locIDstart)
%    [tHit,xHit,xHit_] = extractHits(traj,locIDstart,locIDend)
%
% Inputs:
%    traj - trajectory object
%    locIDstart - (optional) location vector where transition was triggered
%    locIDend - (optional) location vector to where state transitioned
%
% Outputs:
%    tHit - vector of switching times (1 x numSwitches)
%    xHit - array of state vectors before switching (dim_x x numSwitches)
%    xHit_ - array of state vectors after switching (dim_x x numSwitches)
%    idzHit - array of indizes before switching (1 x numSwitches)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: trajectory, hybridAutomaton/simulateRandom

% Authors:       Mark Wetzlinger, Laura Luetzow
% Written:       07-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,3);

% set default values
[locIDstart,locIDend] = setDefaultValues({[],[]},varargin);

% check input arguments
inputArgsCheck({{traj,'att','trajectory'},...
    {locIDstart,'att','numeric'},...
    {locIDend,'att','numeric'}});

% empty object or all continuous-time simulation runs
if all(isempty(traj)) ...
        || all(arrayfun(@(x) isscalar(x.loc) && x.loc == 0,traj,'UniformOutput',true))
    tHit = []; xHit = {}; xHit_ = {}; return
end

% indizes where the locations changes
idx_jump = any(diff(traj.loc,1,2) ~=0, 1);
indizes = 1:length(traj.t)-1;

% shortened code for getting all hits
if isempty(locIDstart) && isempty(locIDend)
    idzHit = indizes(idx_jump);
    tHit = traj.t(idx_jump ~=0);
    xHit = traj.x(:,idx_jump,:);
    xHit_ = traj.x(:,[false idx_jump(1:end-1)],:);
    return
end

% longer version if locations have to be checked
tHit = [];
xHit = [];
xHit_ = [];
idzHit = [];
for j = indizes(idx_jump)

    % check if location before and after jump matches user input
    startcond = isempty(locIDstart) || ...
        (all(size(traj.loc(:,j)) == size(locIDstart)) && ...
        all(traj.loc(:,j) == locIDstart));
    endcond = isempty(locIDend) || ...
        (all(size(traj.loc(:,j+1)) == size(locIDend)) && ...
        all(traj.loc(:,j+1) == locIDend));

    if startcond && endcond
        % time before and after jump are equal so doesn't matter
        tHit  = [tHit traj.t(j)];
        % state vector before jump
        xHit = [xHit traj.x(:,j,:)];
        % state vector after jump
        xHit_ = [xHit_ traj.x(:,j+1,:)];
        % index before jump
        idzHit = [idzHit j];
    end
end

end

% ------------------------------ END OF CODE ------------------------------
