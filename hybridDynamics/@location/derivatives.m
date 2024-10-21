function loc = derivatives(loc,varargin)
% derivatives - compute derivatives for nonlinear reset functions of all or
%    a subset of all transitions within a given location
%
% Syntax:
%    loc = derivatives(loc)
%    loc = derivatives(loc,transIdx)
%    loc = derivatives(loc,transIdx,fpath)
%    loc = derivatives(loc,transIdx,fpath,fname)
%
% Inputs:
%    loc - location object
%    transIdx - index of transition in location object
%    fpath - path to generated file
%    fname - file name
%
% Outputs:
%    loc - updated location object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: transition/derivatives, nonlinearReset/derivatives

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments, set defaults, check values
narginchk(1,4);
[transIdx,fpath,fname] = setDefaultValues({1:numel(loc.transition),...
    [CORAROOT filesep 'models' filesep 'auxiliary' filesep 'location'],...
    'reset'},varargin);
inputArgsCheck({{loc,'att','location','scalar'};...
                {transIdx,'att','numeric','vector'};...
                {fpath,'att','char'};...
                {fname,'att','char'}});

% update transition
for i=1:numel(loc.transition)
    loc.transition(i) = derivatives(loc.transition(i),...
        fpath,sprintf('transition_%i_%s',i,fname));
end

% ------------------------------ END OF CODE ------------------------------
