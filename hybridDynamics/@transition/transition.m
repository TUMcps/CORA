classdef transition
% transition - constructor of the class transition
%
% Syntax:
%    trans = transition()
%    trans = transition(guard,reset,target)
%    trans = transition(guard,reset,target,syncLabel)
%
% Inputs:
%    guard - guard set (contSet object)
%    reset - linearReset or nonlinearReset object
%    target - number of target location
%    syncLabel - synchronization label (parallel hybrid automata only)
%
% Outputs:
%    trans - generated transition object
%
% Example:
%    guard = polytope([0,1],0,[1,0],0);
%    reset = linearReset([0,0;0,0.2]);
%    target = 2;
%    syncLabel = 'on';
%
%    trans = transition(guard,reset,target,syncLabel);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location, hybridAutomaton, parallelHybridAutomaton,
%    linearReset, nonlinearReset

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   30-July-2016
%                10-December-2021 (NK, enable nonlinear reset functions)
%                04-April-2022 (MP, add fields .hasInput/.inputDim to reset)
%                16-June-2022 (MW, add checks for object properties, update handling of reset struct fields)
%                25-June-2024 (TL, added 'real' to symbolic variables)
%                10-October-2024 (MW, update with reset classes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)

    guard;                      % guard set (empty, linear, or nonlinear)
    reset;                      % reset function (linearReset or nonlinearReset)
    target;                     % target location (positive integer)
    syncLabel = '';             % synchronization label (char array)

end

methods
    
    % class constructor
    function trans = transition(varargin)

        % 0. empty
        assertNarginConstructor([0,1,3,4],nargin);
        if nargin == 0
            return
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'transition')
            trans = varargin{1}; return
        end
        
        % 2. parse input arguments: varargin -> vars
        [guard,reset,target,syncLabel] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(guard,reset,target,syncLabel,nargin);

        % 4. compute dependent properties
        [guard,reset,target,syncLabel] = ...
            aux_computeProperties(guard,reset,target,syncLabel);
        
        % 5. assign properties
        trans.guard = guard;
        trans.reset = reset;
        trans.target = target;
        trans.syncLabel = syncLabel;

    end
end

end


% Auxiliary functions -----------------------------------------------------

function [guard,reset,target,syncLabel] = aux_parseInputArgs(varargin)

    % default properties
    guard = []; reset = []; target = []; syncLabel = '';

    % parse arguments
    [guard,reset,target,syncLabel] = setDefaultValues(...
        {guard,reset,target,syncLabel},varargin);

end

function aux_checkInputArgs(guard,reset,target,syncLabel,n_in)

if CHECKS_ENABLED && n_in > 0

    inputArgsCheck({{guard,'att',{'interval','polytope','levelSet','fullspace'}};...
                    {reset,'att',{'abstractReset','struct'},'scalar'};...
                    {target,'att','numeric',{'column','nonempty','integer','positive'}};...
                    {syncLabel,'att','char'}});

    % pre-state dimension of reset function must match dimension of guard
    if isa(reset,'abstractReset') && dim(guard) ~= reset.preStateDim
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Dimension of guard set must match pre-state dimension of reset function.'));
    end

end

end

function [guard,reset,target,syncLabel] = aux_computeProperties(guard,reset,target,syncLabel)

% backward compatibility for structs
if isa(reset,'struct')
    CORAwarning("CORA:deprecated","property","reset","CORA v2025",...
        "Please use a linearReset or nonlinearReset object instead of a struct.",...
        "This change was made to improve code reliability.");

    % try to convert the given struct to a linearReset/nonlinearReset object
    try
        if isfield(reset,'A')
            % linear reset function, has .A and .c, potentially also .B
            if isfield(reset,'B')
                reset = linearReset(reset.A,reset.B,reset.c);
            else
                reset = linearReset(reset.A,[],reset.c);
            end
        elseif isfield(reset,'f')
            % nonlinear reset function
            reset = nonlinearReset(reset.f);
        else
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Reset function must be a linearReset or nonlinearReset object.'));
        end
    catch ME
        % conversion not successful... immediately direct away form struct
        % usage and toward classes
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Reset function must be a linearReset or nonlinearReset object.'));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
