classdef (InferiorClasses = {?interval ?intervalMatrix, ?matZonotope}) mptPolytope < contSet
% mptPolytope - object and copy constructor 
%
% Syntax:  
%    obj = mptPolytope(V)
%    obj = mptPolytope(C,d)
%
% Inputs:
%    V - matrix storing the vertices V = [v_1,...,v_p]^T
%    C - matrix for the inequality constraint C*x <= d
%    d - vector for the inequality constraint C*x <= d
%
% Outputs:
%    obj - generated object
%
% Example: 
%    C = [1 0 -1 0 1; 0 1 0 -1 1]';
%    d = [3; 2; 3; 2; 1];
%    
%    poly = mptPolytope(C,d);
%
%    plot(poly);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    P = [];
    halfspace = [];
end

methods

    function obj = mptPolytope(varargin)

        % If no argument is passed (default constructor)
        if nargin == 0
            disp('mptPolytope needs more input values');
            obj.P=[];
            obj.halfspace=[];

        % If 1 argument is passed
        elseif nargin == 1
            % copy constructor
            if isa(varargin{1},'mptPolytope')
                obj = varargin{1};
            else
                %List elements of the class
                try %MPT V3
                    obj.P=Polyhedron(varargin{1});
                catch %MPT V2
                    obj.P=polytope(varargin{1});
                end
                obj.halfspace=[];
            end

        % If 2 arguments are passed
        elseif nargin == 2
            %List elements of the class
            try %MPT V3
                obj.P=Polyhedron(varargin{1},varargin{2});
            catch %MPT V2
                obj.P=polytope(varargin{1},varargin{2});
            end
            obj.halfspace=[];

        % Else if the parameter is an identical object, copy object    
        elseif isa(varargin{1}, 'mptPolytope')
            obj = varargin{1};

        % Else if not enough or too many inputs are passed    
        else
            disp('This class needs more/less input values');
            obj=[];
        end
        
        % set parent object properties
        obj.dimension = obj.P.Dim;

    end
end

methods (Static = true)
    poly = generateRandom(varargin) % generate random polyope
    poly = enclosePoints(points,varargin) % enclose point cloud with polytope
end

end

%------------- END OF CODE --------------