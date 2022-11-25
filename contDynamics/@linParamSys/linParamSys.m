classdef linParamSys < contDynamics
% linParamSys class (linear parametric system)
%
% Syntax:  
%    obj = linParamSys(A,B)
%    obj = linParamSys(A,B,type)
%    obj = linParamSys(name,A,B)
%    obj = linParamSys(name,A,B,type)
%
% Inputs:
%    A - system matrix
%    B - input matrix
%    name - name of the system
%    type - 'constParam' (constant parameter, default) or 'varParam' (time
%            varying parameter)
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    Ac = [-2 0; 1.5 -3];
%    Aw = [0 0; 0.5 0];
%    A = intervalMatrix(Ac,Aw);
%
%    B = [1; 1];
%
%    sys = linParamSys(A,B,'varParam')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      23-September-2010
% Last update:  01-November-2017 (add possibility to change between constant and varying parameters)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    A = 1;
    B = 0;
    %stepSize = 0.01/max(abs(eig(A.center)));
    constParam = 1; %flag if parameters are constant or time-varying
    stepSize = 1;
    taylorTerms = [];
    mappingMatrixSet = [];
    power = [];
    E = [];
    F = [];
    inputF = [];
    inputCorr = [];
    Rinput = [];
    Rtrans = [];
    RV = [];
    sampleMatrix = [];
end
    
methods
    
    % class constructor
    function obj = linParamSys(varargin)
        
        % default values
        name = 'linParamSys';
        type = 'constParam';
        
        % parse input arguments
        if ischar(varargin{1})
            name = varargin{1};
            A = varargin{2};
            B = varargin{3};
            if nargin > 3
               type = varargin{4}; 
            end
        else
            A = varargin{1};
            B = varargin{2};
            if nargin > 2
                type = varargin{3};
            end
        end
        
        % check input arguments
        if ~ischar(type) || ~ismember(type,{'constParam','varParam'})
           error('Wrong value for input argument "type"!'); 
        end
        
        % number of states and inputs
        states = A.dim;
        if ~isnumeric(B)
            temp = randomSampling(B,1);
            inputs = size(temp,2);
        else
            inputs = size(B,2);
        end
        
         % instantiate parent class
        obj@contDynamics(name,states,inputs,1);
        
        % assign object properties
        obj.A = A;
        obj.B = B;
        
        if strcmp(type,'varParam')
           obj.constParam = 0; 
        end
    end
end
end

%------------- END OF CODE --------------