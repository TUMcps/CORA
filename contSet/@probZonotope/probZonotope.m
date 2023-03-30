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

% Author:       Matthias Althoff
% Written:      03-August-2007 
% Last update:  26-February-2008
%               20-March-2015
%               04-May-2020 (MW, transition to classdef)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % Z ... zonotope matrix
    Z;
    % g ... probabilistic generators
    g;
    % cov ... covariance matrix 
    cov;
    % gauss ... determining if Obj.cov is updated
    gauss;
    % gamma ... cut-off mSigma value
    gamma;
end

methods

    function obj = probZonotope(varargin)

        % parse input
        if nargin > 3
            throw(CORAerror('CORA:tooManyInputArgs',3));
        end
        
        if nargin == 1 && isa(varargin{1}, 'probZonotope')
            % copy constructor  
            obj = varargin{1};
            return
        end

        [Z, g, gamma] = setDefaultValues({[], [], 2}, varargin);
        inputArgsCheck({ ...
            {Z, 'att', 'numeric', 'finite'}; ...
            {g, 'att', 'numeric', 'finite'}; ...
            {gamma, 'att', 'numeric'}; ...
        })

        % assign properties
        obj.Z = Z;
        obj.g = g;
        obj.gamma = gamma;
        
        obj.cov = [];
        obj.gauss = false;
        
        if nargin >= 2
            % update covariance matrix
            obj.cov = sigma(obj);
            obj.gauss = true;
        end
        
    end
end


end

%------------- END OF CODE --------------