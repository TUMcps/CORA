function Zred = reduceCombined(Zbundle,option,varargin)
% reduceCombined - Reduces the order of a zonotope bundle by not reducing
% each zonotope separately, but in a combined fashion
%
% Syntax:  
%    [Zbundle]=reduceCombined(Zbundle,option,varargin)
%
% Inputs:
%    Zbundle - zonotope bundle
%    option - reduction method selector
%    varargin - order and/or filterLength
%
% Outputs:
%    Z - reduced zonotope
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: ---
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      21-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%2 inputs
if nargin==2
    order=1;
    filterLength=[];
%3 inputs
elseif nargin==3
    order=varargin{1};
    filterLength=[];
%4 inputs
elseif nargin==4
    order=varargin{1};
    filterLength=varargin{2};
end


%option='methC'
if strcmp(option,'methC')
    [Zred,t]=reduceMethC(Zbundle,filterLength);   
    

%wrong argument
else
    disp('Error: Second argument is wrong');
    Zred=[];
end


%------------- END OF CODE --------------