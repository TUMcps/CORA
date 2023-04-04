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
%    name - name of the system
%    A - system matrix
%    B - input matrix
%    type - constant/time-varying parameters
%           - 'constParam' (constant parameters, default)
%           - 'varParam' (time-varying parameters)
%
% Outputs:
%    obj - generated linParamSys object
%
% Example:
%    Ac = [-2 0; 1.5 -3];
%    Aw = [0 0; 0.5 0];
%    A = intervalMatrix(Ac,Aw);
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
% Last update:  01-November-2017 (add possibility to change between
%                                 constant and varying parameters)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)

    A = 1;                          % state matrix
    B = 0;                          % input matrix
    constParam logical = true;      % constant or time-varying parameters
    
    % ...
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
            throw(CORAerror('CORA:wrongInputInConstructor',...
                '"Type" has to be either "constParam" or "varParam".'));
        end
        
        % number of states and inputs
        sizeA = dim(A);
        states = sizeA(1);
        if isa(B,'matZonotope') || isa(B,'intervalMatrix')
            sizeB = dim(B);
            inputs = sizeB(2);
        else
            inputs = size(B,2);
        end
        
         % instantiate parent class
        obj@contDynamics(name,states,inputs,1);
        
        % assign object properties
        obj.A = A;
        obj.B = B;
        
        if strcmp(type,'varParam')
            obj.constParam = false; 
        end
    end
end
end

%------------- END OF CODE --------------