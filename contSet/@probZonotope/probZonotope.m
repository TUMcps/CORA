classdef (InferiorClasses = {?interval, ?zonotope}) probZonotope < contSet
% probZonotope - class for probabilistic zonotopes
%
% Syntax:
%    obj = probZonotope(Z,G)
%    obj = probZonotope(Z,G,gamma)
%
% Inputs:
%    Z - zonotope matrix Z = [c,g1,...,gp]
%    G - matrix storing the probabilistic generators G = [g1_, ..., gp_]
%    gamma - cut-off value for plotting. The set is cut-off at 2*sigma,
%            where sigma is the variance
%
% Outputs:
%    obj - Generated object
%
% Example:
%    Z = [10 1 -2; 0 1 1];
%    G = [0.6 1.2; 0.6 -1.2];
%    probZ1 = probZonotope(Z,G);
%    gamma = 1.5;
%    probZ2 = probZonotope(Z,G,gamma);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Authors:       Matthias Althoff
% Written:       03-August-2007 
% Last update:   26-February-2008
%                20-March-2015
%                04-May-2020 (MW, transition to classdef)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    Z;      % zonotope matrix
    g;      % probabilistic generators
    cov;    % covariance matrix

    % internally-set properties
    gauss;  % determining if cov is updated
    gamma;  % cut-off mSigma value
end

methods

    function obj = probZonotope(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'probZonotope')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [Z,g,gamma] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(Z,g,gamma,nargin);

        % 4. assign properties
        obj.Z = Z;
        obj.g = g;
        obj.gamma = gamma;
        obj.gauss = false;

        % 5. compute cov
        obj.cov = sigma(obj);
        obj.gauss = true;
        
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random probZonotope
end

end


% Auxiliary functions -----------------------------------------------------

function [Z,g,gamma] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 3
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % no input arguments
    if nargin == 0
        Z = []; g = []; gamma = [];
        return
    end

    % set default values
    [Z,g,gamma] = setDefaultValues({[],[],2},varargin);

end

function aux_checkInputArgs(Z,g,gamma,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {Z, 'att', 'numeric', 'finite'}; ...
            {g, 'att', 'numeric', 'finite'}; ...
            {gamma, 'att', 'numeric'}; ...
        })
        
    end

end

% ------------------------------ END OF CODE ------------------------------
