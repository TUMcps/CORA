function [E] = plus(obj,S,varargin)
% plus - Overloaded '+' operator for approximating the Minkowski sum of two
% ellipsoids 
%
% Syntax:  
%    E = plus(obj,S)
%    E = plus(obj,S,mode)
%    E = plus(obj,S,L)
%    E = plus(obj,S,L,mode)
%
% Inputs:
%    obj - ellipsoid object
%    S  - set representation (or cell array thereof)
%    L (opt)-directions to use for approximation
%    mode(opt)- type of approximation ('i': inner; 'o': outer)
%
% Outputs:
%    E - ellipsoid after Minkowski sum
%
% Example: 
%
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      09-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
%% parsing & checking
% check if first or second argument is ellipsoid
if ~isa(obj,'ellipsoid') 
    tmp = obj;
    obj = S;
    S = tmp;
end
% parse arguments
if isempty(varargin)
    L = zeros(length(obj.q),0);
    mode = 'o';
elseif length(varargin{1})==1
    if isa(varargin{1},'char')
        mode = varargin{1};
        L = zeros(length(obj.q),0);
    elseif isa(varargin{1},'double')
        L = varargin{1};
        mode = 'o';
    else
        error('Wrong type for 3. input argument!');
    end
elseif length(varargin{1})==2
    if ~isa(varargin{1},'double')
        error('Wrong type for 3. input argument!');
    end
    if ~isa(varargin{1},'char')
        error('Wrong type for 4. input argument!');
    end
    L = varargin{1};
    mode = varargin{2};
else
    error('Wrong number of input arguments!');
end

% handle empty cells, S not a cell etc
S = prepareSetCellArray(S,obj);
if isempty(S)
    E = obj;
    return;
end

% check arguments
if ~any(mode==['i','o'])
   error('mode has to be either "i" (inner approx) or "o" (outer approx)');
end

if size(L,1)~=length(obj.q)
    error('Dimensions of L do not match first argument!');
end

N = length(S);

%% different Minkowski additions
if isa(S{1},'double')
    s = sum(cell2mat(reshape(S,[1,numel(S)])),2);
    E = ellipsoid(obj.Q,obj.q+s);
    return;
end

if isa(S{1},'conPolyZono')
    E = S{1} + obj; 
    for i=2:N
        E = S{i} + E; 
    end
    return; 
end

if isa(S{1},'ellipsoid')
   E = plusEllipsoid(obj,S,L,mode);
   return;
end

error('Set representation is not implemented');
%------------- END OF CODE --------------