classdef ellipsoid < contSet
% ellipsoid - object constructor for ellipsoids
% Some ideas for this class have been extracted from [1].
%
% Description:
%    This class represents ellipsoid objects defined as
%    {x | (x - q)' * Q^(-1) (x - q) <= 1}.
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
%    Obj - generated ellipsoid object
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
% Author:       Victor Gassmann, Matthias Althoff
% Written:      13-March-2019
% Last update:  16-October-2019
%               02-May-2020 (MW, add property validation)
%               29-Mar-2021 (MA, faster eigenvalue computation)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    Q (:,:) {mustBeNumeric,mustBeFinite} = [];
    q (:,1) {mustBeNumeric,mustBeFinite} = [];
    TOL (1,1) double {mustBeNonnegative} = 1e-6;
    dim (1,1) {mustBeInteger} = 0;
    rank (1,1) {mustBeInteger} = 0;
    isdegenerate (1,1) {mustBeNumericOrLogical} = false;
end
   
methods

    function Obj = ellipsoid(varargin)
        
        % (default constructor)
        if nargin==0
            Obj.rank = 0;
            Obj.isdegenerate = true;
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
                %mev = eigs(varargin{1},1,'smallestreal'); % throws warning
                %for ill-conditioned matrix as gives NaN
                mev = min(eig(varargin{1}));
                if ~issymmetric(varargin{1},Obj.TOL) || mev<-Obj.TOL
                    [id,msg] = errConstructor(...
                        'Q needs to be positive semidefinite/symmetric.');
                    error(id,msg);
                end
                Obj.Q = varargin{1};
                if nargin==1
                    Obj.q = zeros(size(Obj.Q,1),1);
                else
                    if length(Obj.Q)~=length(varargin{2})
                        [id,msg] = errConstructor(...
                            'Q and q dimensions are not matching.');
                        error(id,msg);
                    end
                    Obj.q=varargin{2};
                end
                
                %generate parent object
                Obj.dim = dim(Obj);
                Obj.rank = rank(Obj);
                if Obj.rank<Obj.dim
                    Obj.isdegenerate = true;
                end
                
            % Else: too many inputs are passed    
            else
                [id,msg] = errConstructor(...
                    'This class only takes 3 arguments.');
                error(id,msg);
            end 
        end
        
        % set parent object properties
        Obj.dimension = Obj.dim;
        
    end
end
methods (Static=true)
    [E] = generateRandom(varargin)
    E = enclosePoints(points,method) % enclose point cloud with ellipsoid
end
end
%------------- END OF CODE --------------