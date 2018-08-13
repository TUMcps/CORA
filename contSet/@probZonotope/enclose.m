function [pZ] = enclose(varargin)
% enclose - Generates a probabilistic zonotope that encloses two 
% probabilistic zonotopes pZ, A*PZ 
%
% Syntax:  
%    [pZ] = enclose(pZ,A)
%
% Inputs:
%    pZ - first probabilistic zonotope object
%    Ar - system matrix multiplied wit time increment r
%
% Outputs:
%    pZ - probabilistic zonotope, that encloses pZ and A*pZ
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      06-September-2007
% Last update:  03-September-2009
%               04-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

% If 2 arguments are passed
if nargin==2
    pZ1=varargin{1};
    Ar=varargin{2};
    Rtrans=zeros(length(Ar),1);
% If 3 arguments are passed    
elseif nargin==3
    pZ1=varargin{1};
    Ar=varargin{2};
    Rtrans=varargin{3};    
end

%compute enclosure of the probabilistic part
[pZ] = probEnclose(pZ1,Ar,Rtrans);

%------------- END OF CODE --------------