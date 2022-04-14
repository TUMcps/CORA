classdef simResult
% simResult - class that stores simulation results
%
% Syntax:  
%    obj = simResult(x,t)
%    obj = simResult(x,t,loc)
%
% Inputs:
%    x - cell-array storing the simulated trajectories, where each
%        trajectory is a matrix of dimension [N,n]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [N,n]
%    loc - cell-array storing the locations for the simulated trajectories
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet, simulateRandom, simulateRRT

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    x (:,1) {cell} = {};     % states of the simulated trajectories
    t (:,1) {cell} = {};     % time of the simulated trajectories
    loc (:,1) {cell} = {};   % index of the locations (hybrid system only)
end
    
methods
    
    % class constructor
    function obj = simResult(varargin)
        
        % parse input arguments
        if nargin == 2
            obj.x = varargin{1};
            obj.t = varargin{2};
        elseif nargin == 3
            obj.x = varargin{1};
            obj.t = varargin{2};
            obj.loc = varargin{3};
        else
           error('Wrong number of input arguments for class "reachSet"!'); 
        end  
    end   
end
end

%------------- END OF CODE --------------