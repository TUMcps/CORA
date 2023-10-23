classdef testCase
% testCase - class that stores a test case for conformance testing and
% reachset-conformant identification (see Def. 6 in [1]; first publication 
% on this topic was [2]).
%
% Syntax:
%    obj = testCase(y,u,x,dt)
%    obj = testCase(y,u,x0,dt)
%    obj = testCase(y,u,x,dt,name)
%    obj = testCase(y,u,x0,dt,name)
%
% Inputs:
%    y      - (a x q) vector of the measured outputs samples
%    u      - (a x p) vector of input samples
%    x      - (a x n) vector of state samples
%    x0     - (n x 1) vector of initial states
%    dt     - sampling time
%
% Outputs:
%    obj - generated testCase object
%
% References:
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 2022
%    [2] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Stefan Liu
% Written:       15-June-2023             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    y = [];                 % measured outputs
    u = [];                 % inputs
    x = [];                 % states
    initialState = [];      % initialState
    sampleTime = [];        % sample time
    name = [];              % name
    model = [];             % handle to a white box model (if available)
end
    
methods
    
    % class constructor
    function obj = testCase(varargin) 
        
        % check number of input arguments
        if nargin < 4
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 6
            throw(CORAerror('CORA:tooManyInputArgs',8));
        end
        
        obj.y = varargin{1}; % y
        obj.u = varargin{2}; % u
        obj.sampleTime = varargin{4}; % dt
        
        % handle the states
        x = varargin{3}; % x
        if size(x,1)==1 % only one sample
            obj.initialState = x'; % x0
            if size(obj.y,1)==1 % only one measurement
                obj.x = obj.initialState; % state is initial state 
            end
        elseif size(x,2)==1 && size(x,1)~=size(obj.y,1) % not the same number of samples as y
            obj.initialState = x;
        else
            obj.x = x;
            obj.initialState = x(1,:)'; % first row is initial state
        end
        
        % name is additionally specified
        if nargin > 4
            obj.name = varargin{5};
        end
        
        % handle to model is additionally specified
        if nargin > 5
            obj.model = varargin{6};
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
