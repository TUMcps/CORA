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
%    ellipsoid:           Empty set
%    ellipsoid(E):        E-ellipsoids object
%    ellipsoid(Q):        Q positive semidefinite, symmetric matrix
%    ellipsoid(Q,q):      Q, q center of ellipsoid
%    ellipsoid(Q,q,TOL):  Q, q, TOL is tolerance for psd check
%
% Inputs:
%    E - ellipsoid object
%    Q - square shape matrix
%    q - center vector
%
% Outputs:
%    obj - generated ellipsoid object
%
% Example:
%    Q = [2.7 -0.2;-0.2 2.4];
%    q = [1;2];
%    E = ellipsoid(Q, q);
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

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       13-March-2019
% Last update:   16-October-2019
%                02-May-2020 (MW, add property validation)
%                29-March-2021 (MA, faster eigenvalue computation)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    Q;      % shape matrix
    q;      % center
    TOL;    % tolerance
end
   
methods

    function obj = ellipsoid(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'ellipsoid')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [Q,q,TOL] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(Q,q,TOL,nargin);

        % 4. assign properties
        obj.Q = Q;
        obj.q = q;
        obj.TOL = TOL;

    end
end

methods (Static=true)
    E = generateRandom(varargin)
    E = enclosePoints(points,method) % enclose point cloud with ellipsoid
    E = array(varargin)
end

end


% Auxiliary functions -----------------------------------------------------

function [Q,q,TOL] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    if nargin > 3
        % too many input arguments
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % no input arguments
    if nargin == 0
        Q = []; q = []; TOL = 1e-6;
        return
    end

    % assign shape matrix
    Q = varargin{1};
    % set default values
    [q,TOL] = setDefaultValues({zeros(size(Q,1),1),1e-6},varargin(2:end));

end

function aux_checkInputArgs(Q,q,TOL,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {Q, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {q, 'att', 'numeric', {'finite', 'column'}}; ...
            {TOL, 'att', 'numeric', {'nonnegative', 'scalar'}}; ...
        })

        % check dimensions
        if length(Q)~=length(q)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Q and q dimensions are not matching.'));
        end
        mev = min(eig(Q));
        if ~isempty(Q) && (~isApproxSymmetric(Q,TOL) || mev<-TOL)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Q needs to be positive semidefinite/symmetric.'));
        end

    end

end

% ------------------------------ END OF CODE ------------------------------
