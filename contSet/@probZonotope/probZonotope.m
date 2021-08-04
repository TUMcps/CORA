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
%    Obj - Generated Object
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
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % Z ... zonotope matrix
    Z (:,:) {mustBeNumeric,mustBeFinite} = [];
    % g ... probabilistic generators
    g (:,:) {mustBeNumeric,mustBeFinite} = [];
    % cov ... covariance matrix 
    cov (:,:) {mustBeNumeric,mustBeFinite} = [];
    % gauss ... determining if Obj.cov is updated
    gauss (1,1) {mustBeNumericOrLogical} = false;
    % gamma ... cut-off mSigma value
    gamma = 2;
end

methods

    function Obj = probZonotope(varargin)
        
        % default constructor
        if nargin == 0

        % copy constructor  
        elseif nargin == 1 && isa(varargin{1}, 'probZonotope')
            Obj = varargin{1};

        % 2 input arguments
        elseif nargin == 2
            
            % list elements of the class
            Obj.Z = varargin{1}; 
            Obj.g = varargin{2}; 
            Obj.cov = []; 
            Obj.gauss = false; 
            Obj.gamma = 2;      % default value

            % update covariance matrix
            Obj.cov = sigma(Obj);
            Obj.gauss = true;
            
            
        % 3 input arguments
        elseif nargin == 3
            
            % list elements of the class
            Obj.Z = varargin{1}; 
            Obj.g = varargin{2}; 
            Obj.cov = []; 
            Obj.gauss = false; 
            Obj.gamma = varargin{3};

            % update covariance matrix
            Obj.cov = sigma(Obj);
            Obj.gauss = true;

        % error if too many inputs are passed    
        else
            error('Wrong syntax! Type "help probZonotope" for more information.');
        end
        
        % set parent object properties
        Obj.dimension = size(Obj.Z,1);
        
    end
end


end

%------------- END OF CODE --------------