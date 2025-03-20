classdef (InferiorClasses = {?mp}) interval < contSet
% interval - object constructor for real-valued intervals 
%
% Description:
%    This class represents interval objects defined as
%    {x | a_i <= x <= b_i, \forall i = 1,...,n}.
%
% Syntax:
%    obj = interval(I)
%    obj = interval(a)
%    obj = interval(a,b)
%
% Inputs:
%    I - interval object
%    a - lower limit
%    b - upper limit
%
% Outputs:
%    obj - generated interval object
%
% Example:
%    a = [1;-1];
%    b = [2;3];
%    I = interval(a,b);
%    plot(I,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       19-June-2015
% Last update:   18-November-2015
%                26-January-2016
%                15-July-2017 (NK)
%                01-May-2020 (MW, delete redundant if-else)
%                20-March-2021 (MW, error messages)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                29-March-2023 (TL, optimized constructor)
%                08-December-2023 (MW, handle [-Inf,-Inf] / [Inf,Inf] case)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    inf;    % lower bound
    sup;    % upper bound
end

methods

    % class constructor
    function obj = interval(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end
        assertNarginConstructor(1:2,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'interval')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [lb,ub] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(lb,ub,nargin);

        % 4. compute properties (deal with corner cases)
        [lb,ub] = aux_computeProperties(lb,ub);

        % 5. assign properties
        obj.inf = lb;
        obj.sup = ub;

        % 6. set precedence (fixed)
        obj.precedence = 120;

    end
    

    function ind = end(obj,k,n)
    % overloads the end operator for referencing elements, e.g. I(end,2),
        ind = size(obj,k);
    end
end

methods (Static = true)
    I = generateRandom(varargin) % generates random interval
    I = enclosePoints(points) % enclosure of point cloud
    I = empty(n) % instantiates an empty interval
    I = Inf(n) % instantiates a fullspace interval
    I = origin(n) % instantiates an interval representing the origin in R^n
end

methods (Access = protected)
    [abbrev,printOrder] = getPrintSetInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [lb,ub] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables
    
    % assign lower and upper bound
    [lb,ub] = setDefaultValues({[],[]},varargin);

    % set upper bound to value of lower bound if only one value given
    if isnumeric(lb) && ~isempty(lb) && isempty(ub)
        ub = lb;
    end

end

function aux_checkInputArgs(lb,ub,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0
        
        inputArgsCheck({ ...
            {lb, 'att', 'numeric', 'nonnan'}; ...
            {ub, 'att', 'numeric', 'nonnan'}; ...
        })

        if ~isempty(lb) && ~isempty(ub)
            if ~all(size(lb) == size(ub))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Limits are of different dimension.'));
            % elseif length(size(lb)) > 2
            %     throw(CORAerror('CORA:wrongInputInConstructor',...
            %         'Only 1d and 2d intervals are supported.'));
            elseif ~all(lb <= ub,"all")
                % check again using tolerance (little bit slower)
                % tolerance is chosen as 1e-6 since many interval
                % enclosures are computed using linear programs
                if ~all(lb < ub | withinTol(double(lb),double(ub),1e-6),"all")
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Lower limit larger than upper limit.'));
                end
            end
        end
        
    end

end

function [lb,ub] = aux_computeProperties(lb,ub)
% if one dimension is [-Inf,-Inf] or [Inf,Inf], the interval is empty
% (cannot be displayed for interval matrices -> throws error)

    % assign correct size if empty interval
    if isempty(lb) && isempty(ub)
        ub = zeros(size(lb));
    end

    % check if given interval is empty
    if any(any( isinf(lb) & isinf(ub) & (sign(lb) == sign(ub)) ))
        n = size(lb);
        if all(n > 1)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Empty interval matrix cannot be instantiated'));
        end
        n(n==1) = 0;
        lb = zeros(n);
        ub = zeros(n);
    end

end

% ------------------------------ END OF CODE ------------------------------
