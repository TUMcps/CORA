function [E] = minus(obj,S,varargin)
% plus - Overloaded '-' operator for approximating the Minkowski difference
% of ellipsoids
%
% Syntax:  
%    E = minus(obj,S)
%    E = minus(obj,S,mode)
%    E = minus(obj,S,L)
%    E = minus(obj,S,L,mode)
%
% Inputs:
%    obj - ellipsoid object
%    S  - set representation (or cell array thereof)
%    L (opt)-directions to use for approximation
%    mode(opt)- type of approximation ('i': inner; 'o': outer)
%
% Outputs:
%    E - ellipsoid after Minkowski difference
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
    if ~isa(varargin{2},'char')
        error('Wrong type for 4. input argument!');
    end
    L = varargin{1};
    mode = varargin{2};
else
    error('Wrong number of input arguments!');
end

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
    error('Dimension of L does not match dimension of ellipsoid!');
end

if isa(S,'double')
    if any(size(S)==[center(obj),1])
        error('Wrong size for second argument of type "double"');
    end
end

N = length(S);
%% different Minkowski differences

if isa(S{1},'double')
    s = sum(cell2mat(reshape(S,[1,numel(S)])),2);
    E = ellipsoid(obj.Q,obj.q-s);
    return;
end

if isa(S{1},'ellipsoid')
    E = minusEllipsoid(obj,S{1},L,mode);
    for i=2:N
        E = minusEllipsoid(E,S{i},L,mode);
    end
    return;
end
%------------- END OF CODE --------------