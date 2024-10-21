function HA = derivatives(HA,varargin)
% derivatives - compute derivatives for nonlinear reset functions of all
%    transitions of all or a subset of all locations within a hybrid
%    automaton
%
% Syntax:
%    HA = derivatives(HA)
%    HA = derivatives(HA,locIdx)
%    HA = derivatives(HA,locIdx,fpath)
%    HA = derivatives(HA,locIdx,fpath,fname)
%
% Inputs:
%    HA - hybridAutomaton object
%    locIdx - indices of locations
%    fpath - path to generated file
%    fname - file name
%
% Outputs:
%    HA - updated hybridAutomaton object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/derivatives, transition/derivatives,
%    nonlinearReset/derivatives

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments, set defaults, check values
narginchk(1,4);
[locIdx,fpath,fname] = setDefaultValues({1:numel(HA.location),...
    [CORAROOT filesep 'models' filesep 'auxiliary' filesep ...
    'hybridAutomaton' filesep HA.name],...
    'nonlinear_reset_function'},varargin);
inputArgsCheck({{HA,'att','hybridAutomaton','scalar'};...
                {locIdx,'att','numeric','vector'};...
                {fpath,'att','char'};...
                {fname,'att','char'}});

% loop over all locations
for i=1:numel(HA.location)
    % make a new folder for each location (note: we do not use the name of
    % the location since there may be duplicates)
    HA.location(i) = derivatives(HA.location(i),...
        numel(HA.location(i).transition),...
        sprintf('%s%slocation_%i',fpath,filesep,i),fname);
end

% ------------------------------ END OF CODE ------------------------------
