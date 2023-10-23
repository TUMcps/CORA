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
%    reset - reset function (struct)
%            - linear (without inputs): fields 'A','c'
%                x_ = Ax + c, where x_ is the state after reset
%            - linear (with inputs): fields 'A','B','c'
%                x_ = Ax + Bu + c, where x_ is the state after reset
%            - nonlinear (with/without inputs): field 'f' (function handle)
%    target - number of target location
%    syncLabel - synchronization label (only for use in parallel hybrid
%                automata)
%
% Outputs:
%    trans - generated transition object
%
% Example:
%    guard = conHyperplane([1,0],0,[0,1],0);
%    reset = struct('A',[0, 0; 0, 0.2],'c',zeros(2,1));
%    target = 2;
%    syncLabel = 'on';
%
%    trans = transition(guard,reset,target,syncLabel);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, location, parallelHybridAutomaton

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   30-July-2016
%                10-December-2021 (NK, enable nonlinear reset functions)
%                04-April-2022 (MP, add fields .hasInput/.inputDim to reset)
%                16-June-2022 (MW, add checks for object properties, update handling of reset struct fields)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % guard set (empty, linear, or nonlinear)
    guard = [];

    % reset function
    % - linear: struct with fields 'A'(,'B'),'c' for Ax (+ Bu) + c
    % - nonlinear: struct with field 'f'
    reset struct = struct();

    % target location
    target {mustBeNumeric} = [];

    % synchronization label
    syncLabel (1,:) char = '';
end

methods
    
    % class constructor
    function trans = transition(varargin)

        % parse input arguments
        if nargin == 0
            return
        elseif nargin == 1 && isa(varargin{1},'transition')
            % copy constructor
            trans = varargin{1}; return
        elseif nargin < 3
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
            
        % assign object properties
        trans.guard = varargin{1};
        trans.reset = varargin{2};
        trans.target = varargin{3};
        if nargin == 4
            % 4. check synchronization label already here (after assignment
            % other types are also accepted as char for some reason)
            if ~ischar(varargin{4})
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Synchronization label has to be a char-array.'));
            end
            trans.syncLabel = varargin{4};
        end
        
        % 1. check guard set: either empty or admissible set representation
        guardlist = {'interval','polytope','levelSet','conHyperplane','fullspace'};
        if ~any(ismember(class(trans.guard),guardlist))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Property ''guard'' has to be one of the following classes:\n'...
                '  ' strjoin(guardlist,', ') '.']));
        end

        % 2. check reset function: correct fields, add internal fields
        % 'hasInput' and 'inputDim', compute derivatives for nonlinear f
        list = fields(trans.reset);

        % check fields of reset struct:
        % - from user instantiation: A,B,c,f
        % - internal instantations: hasInput, inputDim, stateDim, J, Q, T
        admissibleFields = {'A','B','c','f',...
            'hasInput','inputDim','stateDim','J','Q','T'};
        givenFields = ismember(admissibleFields,list);
        redundantFields = rmfield(trans.reset,admissibleFields(givenFields));
        if ~isempty(fields(redundantFields))
            % throw error if any fields are redundant
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['transition.reset has redundant fields: '...
                strjoin(fields(redundantFields),', '),'.']));
        end
        
        % assign reset function properties
        if ~all(isfield(trans.reset,{'stateDim','inputDim','hasInput'}))
            % assumption: reset function has input arguments @(x,u)
            trans.reset.hasInput = false;
            trans.reset.inputDim = 0;
            if isfield(trans.reset,'f')
                % nonlinear reset function x_ = f(x,u)
                % read out length of input arguments to evaluate reset function
                [temp,trans.reset.stateDim] = inputArgsLength(trans.reset.f,2);
                % size of inputs to reset function
                trans.reset.inputDim = temp(2);
                trans.reset.hasInput = trans.reset.inputDim > 0;
            elseif isfield(trans.reset,'B')
                % linear reset function with inputs: x_ = Ax + Bu + c
                trans.reset.hasInput = true;
                trans.reset.stateDim = size(trans.reset.B,1);
                trans.reset.inputDim = size(trans.reset.B,2);
            else
                % linear reset function without inputs
                trans.reset.stateDim = size(trans.reset.A,1);
            end
        end

        % check whether matrices of linear transitions have correct sizes
        aux_checkLinearReset(trans);
        
        % pre-compute derivatives for nonlinear resets: stored in fields
        % .J (Jacobian), .Q (Hessian), .T (third-order tensor)
        if isfield(trans.reset,'f') && ~all(isfield(trans.reset,{'J','Q','T'}))
            % skip computation if J, Q, and T have already been computed
            % (this occurs during internal re-instantation of transitions)
            trans = aux_compDerivatives(trans);
        end

        % 3. check target: non-empty, integer, larger than zero
        if isempty(trans.target) || any(mod(trans.target,1) ~= 0) || any(trans.target < 1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'All targets of a transition have to be integer values greater than zero.'));
        end

    end
    

    % Auxiliary functions -------------------------------------------------

    function aux_checkLinearReset(trans)
    % ensure that linear reset function x_ = Ax + Bu + c is properly
    % defined

        if CHECKS_ENABLED && isfield(trans.reset,'A')
            % size of next state vector and current state vector
            if isfield(trans.reset,'c')
                if ~isvector(trans.reset.c)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Constant offset in linear reset function has to be a vector.'));
                end

                if size(trans.reset.A,1) ~= length(trans.reset.c)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Linear reset function: Mismatch in rows of A and length of c.'));
                end
            end

            if trans.reset.hasInput && size(trans.reset.A,1) ~= size(trans.reset.B,1)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Linear reset function: Mismatch in rows of A and rows of B.'));
            end
        end

    end

    function trans = aux_compDerivatives(trans)
    % pre-compute all derivatives required for the enclosure of the 
    % nonlinear reset function by a Taylor series expansion of order 2
    % additional fields: .J (Jacobian), .Q (Hessian), .T (third-order tensor)

        % state dimension and input dimension of reset function
        n = trans.reset.stateDim;
        m = trans.reset.inputDim;

        % create symbolic variables for reset function
        z = sym('x',[n,1]);
        u = sym('u',[m,1]);

        % if reset function depends on inputs, we create the substitute
        % state z = [x;u] which is used for the rest of this function
        if trans.reset.hasInput
            % extended state z = [x;u]
            z = [z;u];
            % alter reset function to accept single input z
            trans.reset.f = @(z) trans.reset.f(z(1:n),z(n+1:n+m));
        end
        
        % evaluate f symbolically
        f = trans.reset.f(z);
        
        % first-order derivative
        J = jacobian(f,z);
        trans.reset.J = matlabFunction(J,'Vars',{z});
        
        % check if Jacobian is constant:
        % insert NaN for all variables -> all entries with variables become
        % NaN, all others simple double
        Jconst = ~isnan(double(subs(J,z,NaN(length(z),1))));
        
        % second-order derivative
        trans.reset.Q = cell(n,1);
        for i = 1:n
            if all(Jconst(i,:))
                % skip linear Jacobians
                trans.reset.Q{i,1} = 0;
            else
                % compute Hessian of f
                trans.reset.Q{i,1} = matlabFunction(hessian(f(i),z),'Vars',{z});
            end
        end
        
        % third-order derivative
        trans.reset.T = cell(n,n+m);
        for i = 1:n
            for j = 1:n+m
                if Jconst(i,j)
                    % skip linear Jacobians
                    trans.reset.T{i,j} = 0;
                else
                    % third-order tensor = Hessian of Jacobian
                    trans.reset.T{i,j} = matlabFunction(hessian(J(i,j),z),'Vars',{z});
                end
            end
        end
    end

end
end

% ------------------------------ END OF CODE ------------------------------
