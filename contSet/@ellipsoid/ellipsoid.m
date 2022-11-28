classdef ellipsoid < contSet
% ellipsoid - object constructor for ellipsoids
% Some ideas for this class have been extracted from [1].
%
% Description:
%    This class represents ellipsoid objects defined as
%    {x | (x - q)' * Q^(-1) (x - q) <= 1}
%    in the non-degenerate case, and are defined using the support function
%    of ellipsoids (see ellipsoid/supportFunc.m for details) in the general
%    case.
%
% Syntax:  
%    object constructor: Obj = ellipsoids(varargin)
%    copy constructor: Obj = otherObj
%
% (Overloaded) Methods:
%    - ellipsoid:           Empty set
%    - ellipsoid(E):        E-ellipsoids object
%    - ellipsoid(Q):        Q positive semidefinite, symmetric matrix
%    - ellipsoid(Q,q):      Q, q center of ellipsoid
%    - ellipsoid(Q,q,TOL):  Q, q, TOL is tolerance for psd check
%
% Inputs:
%    E - ellipsoid object
%    Q - square shape matrix
%    q - center vector
%
% Outputs:
%    obj - generated ellipsoid object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference on
%       Decision and Control (pp. 1498-1503). IEEE.

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
end
   
methods

    function obj = ellipsoid(varargin)
        
        if nargin >= 1
        % (copy constructor)
            if nargin == 1 && isa(varargin{1},'ellipsoid')
                obj = varargin{1};
            
            % If 1-3 arguments are passed
            elseif nargin >=1 && nargin<=3
                %check TOL first as needed later
                if nargin==3
                    obj.TOL = varargin{3};
                end
                

                obj.Q = varargin{1};
                if nargin==1
                    obj.q = zeros(size(obj.Q,1),1);
                else
                    if length(obj.Q)~=length(varargin{2})
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Q and q dimensions are not matching.'));
                    end
                    obj.q=varargin{2};
                end

                % if Q and q have correct dimensions, but are empty use
                % default constructor
                if isempty(obj.Q) && isempty(obj.q)
                    obj = ellipsoid;
                    return;
                end

                %Check psd and symmetry of Q
                %mev = eigs(varargin{1},1,'smallestreal'); % throws warning
                %for ill-conditioned matrix as gives NaN
                mev = min(eig(obj.Q));
                if ~isApproxSymmetric(obj.Q,obj.TOL) || mev<-obj.TOL
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Q needs to be positive semidefinite/symmetric.'));
                end
                
            % Else: too many inputs are passed    
            else
                throw(CORAerror('CORA:tooManyInputArgs',3));
            end 
        end
        
        % set parent object properties
        obj.dimension = dim(obj);
        
    end
end

methods (Static=true)
    E = generateRandom(varargin)
    E = enclosePoints(points,method) % enclose point cloud with ellipsoid
    E = array(varargin)
end

end
%------------- END OF CODE --------------