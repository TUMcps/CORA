function varargout = conform(varargin)
% conform - performs reachset conformance by computing the reachable 
%    of a system and checking whether all measurements are included.
%    Optionally, conformance checking can be computed by passing the 'check'
%    parameter
%
% Syntax:
%    [params, R, simRes, unifiedOutputs, updatedSys] = conform(sys,'synth',params,options)
%    [res, R, simRes, unifiedOutputs] = conform(sys,'check',params,options)
%
% Inputs:
%    sys - contDynamics object
%    type - 'synth' or 'check', can be omitted for synthesis
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    params - struct with conformant parameters
%    R - reachSet object (only time steps for which measurments exist)
%    simRes - states of the rapidly exploring random tree
%    unifiedOutputs - unified test cases (only deviation to nominal solution)
%    updatedSys - updated system
%    res - result: true/false  (for conformance checking)
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       15-June-2023
% Last update:   12-October-2023 (TL, internal split confSynth/confCheck)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[sys,type,params,options] = aux_parseInput(varargin{:});

% 2. decide type
switch(type)
    case 'synth' % conformance synthesis
        varargout = cell(1, 5); % 5 output arguments
        [varargout{:}] = confSynth(sys,params,options);

    case 'check' % conformance checking
        varargout = cell(1, 4); % 4 output arguments
        [varargout{:}] = confCheck(sys,params,options);

    otherwise
        throw(CORAerror('CORA:wrongValue','second','type needs to be ''synth'' or ''check''.'))
end

end


% Auxiliary functions -----------------------------------------------------

function [sys,type,params,options] = aux_parseInput(varargin)
    % only rough checks here, checks are then done in 
    % confSynth/confCheck
    if nargin < 3
        throw(CORAerror('CORA:notEnoughInputArgs',3))
    elseif nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4))
    end

    % read parameters
    sys = varargin{1};
    if nargin == 4 % conform(sys, type, params, options)
        type = varargin{2};
        params = varargin{3};
        options = varargin{4};

    elseif nargin == 3 % conform(sys, params, options)
        type = 'synth';
        params = varargin{2};
        options = varargin{3};
    end

end

% ------------------------------ END OF CODE ------------------------------
