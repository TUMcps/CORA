function I = minkDiff(I,S,varargin)
% minkDiff - compute the Minkowski difference of two intervals:
%            I1 - I2 = I <-> I + I2 \subseteq I1
%
% Syntax:
%    I = minkDiff(I,S)
%    I = minkDiff(I,S,type)
%
% Inputs:
%    I - interval object
%    S - interval object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'inner')
%
% Outputs:
%    I - interval object after Minkowski difference
%
% Example: 
%    I1 = interval([-2;-1],[3;3]);
%    I2 = interval([-1;-1],[1;1]);
%
%    I = minkDiff(I1,I2);
%
%    figure; hold on;
%    plot(I1);
%    plot(I2,[1,2],'r');
%    plot(I,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Authors:       Niklas Kochdumper
% Written:       10-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
type = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{I,'att','interval'};
                {S,'att',{'interval','contSet','numeric'}};
                {type,'str',{'exact','inner'}}});

% different algorithms for different set representations
if isnumeric(S)
    
   I = I + (-S);
   
elseif isa(S,'interval')
    try
       I = interval(infimum(I)-infimum(S),supremum(I)-supremum(S)); 
    catch
       I = []; 
    end
    
else
    
    % parse input arguments
    type = 'approx';
    if nargin > 2 && ~isempty(varargin{1})
        type = varargin{1};
    end
    
    % check input arguments
    if strcmp(type,'exact')
        throw(CORAerror('CORA:noExactAlg',I,S));
    end
    
    % compute inner-approximation of the Minkowski difference
    infi = infimum(I); sup = supremum(I); n = dim(I);
    
    for i = 1:n
       temp = zeros(n,1);
       temp(i) = 1;
       sup(i) = sup(i) - supportFunc_(S,temp,'upper','interval',8,1e-3);
       infi(i) = infi(i) + supportFunc_(S,-temp,'upper','interval',8,1e-3);
    end
    
    % construct resulting interval
    try
        I = interval(infi,sup);
    catch
        I = []; 
    end
end

% ------------------------------ END OF CODE ------------------------------
