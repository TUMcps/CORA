function [V] = vertices(varargin)
% vertices - Returns potential vertices of a zonotope
% WARNING: Do not use this function for high order zonotopes -
% computational complexity grows exponential!
%
% Syntax:  
%    [V] = vertices(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    V=vertices(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervallhull,  polytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  30-January-2008
%               23-June-2009
%               24-June-2010
%               11-July-2012
%               30-July-2016
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin==1
    Z=varargin{1};
    dim=[];
    option='';
elseif nargin==2
    Z=varargin{1};
    dim=varargin{2};
    option='';
elseif nargin==3
    Z=varargin{1};
    dim=varargin{2};
    option='alternative';
end

%get matrix from object
if ~isempty(dim)
    [Z] = project(Z,dim);
    Z=deleteZeros(Z);
end

%convert to polytope
if ~strcmp(option,'alternative')
    try
        options.polytopeType='';
        P=polytope(Z,options);
        %compute extreme points
        V=vertices(P);
        if isempty(V)
            option='alternative';
        end
    catch
        try
            options.polytopeType='mpt';
            P=polytope(Z,options);
            %compute extreme points
            V=vertices(P);
            if isempty(V)
                option='alternative';
            end
        catch
            option='alternative';
        end
    end
end

if strcmp(option,'alternative')
    vertexArray=alternative(Z.Z);
    %create vertices object
    V=vertices(vertexArray);
end




function vertexArray=alternative(Z);

%first vertex is the center of the zonotope
vertexArray=Z(:,1);

%Generate further potential vertices in the loop
for iVertex=1:length(Z(1,2:end))
    translation=Z(:,iVertex+1)*ones(1,length(vertexArray(1,:)));
    V=[vertexArray+translation,vertexArray-translation];
    %remove inner points
    if iVertex>length(Z(:,1))
        try
            K = convhulln(V');
            indices = unique(K);
            vertexArray=V(:,indices);
        catch
            disp('Convex hull failed')
            vertexArray=V;
        end
    else
        vertexArray=V;
    end
end

%------------- END OF CODE --------------