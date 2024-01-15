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
%    ellipsoid(E)
%    ellipsoid(Q)
%    ellipsoid(Q,q)
%    ellipsoid(Q,q,TOL)
%
% Inputs:
%    E - ellipsoid object
%    Q - square, positive semi-definite shape matrix
%    q - center vector
%    TOL - tolerance
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
% See also: none
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

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'ellipsoid')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [Q,q,TOL] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(Q,q,TOL);

        % 4. compute properties
        [Q,q,TOL] = aux_computeProperties(Q,q,TOL);

        % 5. assign properties
        obj.Q = Q;
        obj.q = q;
        obj.TOL = TOL;

    end
end

methods (Static=true)
    E = generateRandom(varargin) % generates a random ellipsoid
    E = enclosePoints(points,method) % enclose point cloud with ellipsoid
    E = array(varargin)
    E = empty(n) % instantiates an empty ellipsoid
end

end


% Auxiliary functions -----------------------------------------------------

function [Q,q,TOL] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    if nargin > 3
        % too many input arguments
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % assign shape matrix
    Q = varargin{1};
    % set default values
    [q,TOL] = setDefaultValues({zeros(size(Q,1),1),1e-6},varargin(2:end));

end

function aux_checkInputArgs(Q,q,TOL)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED


        % allow empty Q matrix for ellipsoid.empty
        if isempty(Q)
            % only ensure that q is also empty
            if ~isempty(q)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Shape matrix is empty, but center is not.'));
            end
            return
        end

        inputArgsCheck({ ...
            {Q, 'att', 'numeric', {'finite', 'matrix'}}; ...
            {q, 'att', 'numeric', {'finite', 'column'}}; ...
            {TOL, 'att', 'numeric', {'nonnegative', 'scalar'}}; ...
        });

        % shape matrix needs to be square
        if size(Q,1) ~= size(Q,2)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The shape matrix needs to be a square matrix.'));
        end

        % check dimensions
        if ~isempty(q) && length(Q) ~= size(q,1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Dimensions of the shape matrix and center are different.'));
        end
        mev = min(eig(Q));
        if ~isempty(Q) && (~isApproxSymmetric(Q,TOL) || mev<-TOL)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'The shape matrix needs to be positive semidefinite/symmetric.'));
        end

    end

end

function [Q,q,TOL] = aux_computeProperties(Q,q,TOL)

    if isempty(Q)
        Q = zeros(0,0);
        if ~isempty(q)
            q = zeros(0,0);
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
