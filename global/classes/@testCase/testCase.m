classdef testCase
% testCase - class that stores a test case for conformance testing and
%    reachset-conformant identification (see Def. 6 in [1]; first
%    publication on this topic was [2]).
%
% Syntax:
%    obj = testCase(y,u,x,dt)
%    obj = testCase(y,u,x0,dt)
%    obj = testCase(y,u,x,dt,name)
%    obj = testCase(y,u,x0,dt,name)
%
% Inputs:
%    y - (a x q x s) vector of the measured outputs samples
%    u - (a x p x s) vector of input samples
%    x - (a x n) vector of state samples
%    x0 - (n x 1 x s) vector of initial states
%    dt - sampling time
%    name - name of the test case
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
% See also: none

% Authors:       Matthias Althoff, Stefan Liu
% Written:       15-June-2023             
% Last update:   25-July-2023 (LL, add uncertain initial state)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    y = [];                 % measured outputs
    y_a = [];               % deviation of measured outputs from the nominal
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
        assertNarginConstructor(4:6,nargin);
        
        obj.y = varargin{1}; % y

        % transform measurement cell array to 3D array
        if iscell(obj.y)
            y = zeros(size(obj.y{1},1),size(obj.y{1},2),length(obj.y));
            for i=1:length(obj.y)
                y(:,:,i) = obj.y{i};
            end
            obj.y = y;
        end
        obj.u = varargin{2}; % u
        obj.sampleTime = varargin{4}; % dt
        
        % handle the states
        x = varargin{3}; % x
        if size(x,2) == 1 % isa(x,'contSet') % initial state is uncertain
            obj.initialState = x;
        else
            obj.x = x;
            obj.initialState = permute(x(1,:,:), [2 1 3]); % first row is initial state
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


    function obj = set_u(obj,u)
        obj.u = u;
    end

    function obj = compute_ya(obj,sys)
        % compute mesurement deviation y_a by subtracting the nominal 
        % solution from the measurement trajectory

        y_nom = zeros(size(obj.y));
        params.tFinal = sys.dt * size(obj.y,1) - sys.dt;
        for s =1:size(obj.y,3)
            params.x0 = obj.initialState(:,:,s);
            params.u = obj.u(:,:,s)';
            [~,~,~,y_nom(:,:,s)] = simulate(sys, params);
        end
        obj.y_a = obj.y - y_nom;
    end   

    function obj = reduceLength(obj,n_k)
        obj.u = obj.u(1:n_k,:,:);
        obj.y = obj.y(1:n_k,:,:);
    end

    function obj = combineTestCases(obj,obj2)
        % combine two test cases in one test case to simplify computations 
        % for linear systems

        % stack intial states
        obj.initialState = cat(3, obj.initialState, obj2.initialState);    

        % stack input trajectories, fill with NaN if different length  
        u_diff = size(obj.u,1) - size(obj2.u,1);
        dim_u = size(obj.u,2);
        obj.u = cat(3, [obj.u; NaN(max(-u_diff,0),dim_u,size(obj.u,3))], ...
            [obj2.u; NaN(max(u_diff,0),dim_u,size(obj2.u,3))]);

        % stack output trajectories, fill with NaN if different length  
        y_diff = size(obj.y,1) - size(obj2.y,1);
        dim_y = size(obj.y,2);
        obj.y = cat(3, [obj.y; NaN(max(-y_diff,0),dim_y,size(obj.y,3))], ...
            [obj2.y; NaN(max(y_diff,0),dim_y,size(obj2.y,3))]);
    end

    function testSuite = setInitialStateToMeas(obj,p, tol)
        % set initial state to the first p measurements (necessary for
        % input.output models)
        if nargin < 3
            tol = 1e-12;
        end

        if mean(abs(diff(obj.y(1:p,:,:), 1, 3)),'all') < tol
            % all measurement trajectories start with the same measurements
            obj.initialState = reshape((obj.y(1:p,:,1))',[],1);
            testSuite{1} = obj;
        else
            % generate new test cases for each measurements trajectory
            n_s = size(obj.y,3);
            testSuite = cell(n_s,1);
            for s=1:n_s
                testSuite{s} = obj;
                testSuite{s}.initialState = reshape((obj.y(1:p,:,s))',[],1);
                testSuite{s}.y = obj.y(:,:,s);
            end
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
