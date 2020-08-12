classdef ellipsoid
% ellipsoid - Object and Copy Constructor 
% Some ideas for this class have been extracted from [1].
%
% Syntax:  
%    object constructor: Obj = ellipsoids(varargin)
%    copy constructor: Obj = otherObj
%
%(Overloaded) Methods:
%    - ellipsoid:           Empty set
%    - ellipsoid(E):        E-ellipsoids object
%    - ellipsoid(Q):        Q positive semidefinite, symmetric matrix
%    - ellipsoid(Q,q):      Q, q center of ellipsoid
%    - ellipsoid(Q,q,TOL):  Q, q, TOL is tolerance for psd check
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%
% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  16-October-2019
%               02-May-2020 (MW, add property validation)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    Q (:,:) {mustBeNumeric,mustBeFinite} = [];
    q (:,1) {mustBeNumeric,mustBeFinite} = [];
    contSet = [];
    TOL (1,1) double {mustBeNonnegative} = 1e-12;
    dim (1,1) {mustBeInteger} = 0;
    isdegenerate (1,1) {mustBeNumericOrLogical} = false;
end
   
methods

    function Obj = ellipsoid(varargin)
        
        % (default constructor)
        if nargin==0
            Obj.contSet = contSet();
            
        else
        % (copy constructor)
            if nargin == 1 && isa(varargin{1},'ellipsoid')
                Obj = varargin{1};
            
            % If 1-3 arguments are passed
            elseif nargin >=1 && nargin<=3
                %check TOL first as needed later
                if nargin==3
                    Obj.TOL = varargin{3};
                end

                %Check psd and symmetry of Q
                mev = min(eig(varargin{1}));
                if ~issymmetric(varargin{1},Obj.TOL) || mev<-Obj.TOL
                    error('Q needs to be positive semidefinite/symmetric');
                end
                Obj.Q=varargin{1};
                if nargin==1
                    Obj.q = zeros(size(Obj.Q,1),1);
                else
                    if length(Obj.Q)~=length(varargin{2})
                        error('Q and q dimensions are not matching');
                    end
                    Obj.q=varargin{2};
                end
                
                %generate parent object
                Obj.contSet = contSet(length(Obj.Q));
                %Determine the rank of Q; helpful to e.g. determine whether
                %ellipsoid is degenerate or not
                Obj.dim=rank(Obj.Q,Obj.TOL);
                if Obj.dim<length(Obj.Q)
                    Obj.isdegenerate=true;
                end
                
            % Else: too many inputs are passed    
            else
                error('This class only takes 3 arguments');
            end
            
        end
        
    end
end
methods (Static=true)
    [E] = generateRandom(varargin)
    E = enclosePoints(points) % enclose point cloud with ellipsoid
end
end
%------------- END OF CODE --------------