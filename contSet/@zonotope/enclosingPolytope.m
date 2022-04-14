function P = enclosingPolytope(Z,varargin)
% enclosingPolytope - converts a zonotope to a polytope representation
%
% Syntax:  
%    P = enclosingPolytope(Z,varargin)
%
% Inputs:
%    Z - zonotope object
%    options - options containing method of enclosure
%
% Outputs:
%    P - polytope object
%
% Example: 
%    Z = zonotope([0;0],rand(2,5));
%    P = enclosingPolytope(Z);
%    plot(Z); hold on;
%    plot(P,[1,2],'r');
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      18-September-2007
% Last update:  26-August-2010
%               20-October-2010
%               25-July-2016 (intervalhull replaced by interval)
%               05-April-2017
%               20-April-2018 (methC replaced by PCA)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin==1
     %filterLength=[];
     tightConversion=1;
     options.polytopeType = 'mpt';
elseif nargin==2
     options=varargin{1};
     
     
%      %filter length
%      if isfield(options,'filterLength')
%          filterLength=options.filterLength;
%      else
%          dim=length(Z.Z(:,1));
%          filterLength(1)=dim+2;
%          filterLength(2)=dim+3;
%      end
     
     %flag for tight conversion
     if isfield(options,'tightPolytopeConversion')
         tightConversion = options.tightPolytopeConversion;
     else 
         tightConversion = true;
     end
end


% %use simple reduction
% if isempty(filterLength)
%      Zred=zonotope(interval(Z));
%      P=parallelotope(Zred, options);
%      
% %filter length is specified     
% else
     if tightConversion
         %solution1 (axis-aligned):
         Zred=zonotope(interval(Z));
         P=polytope(Zred, options); 
         %solution 2 (method C):
         %Zred=reduce(Z,'methC',1,filterLength);
         Zred=reduce(Z,'pca');
         Zred = repair(Zred,Z);
         Padd = polytope(Zred, options);
         %intersect results
         P=P&Padd;
     else
         %solution 1 (method C):
         %Zred1 = reduce(Z,'methC',1,filterLength);
         Zred1 = reduce(Z,'pca');
         Zred1 = repair(Zred1,Z);
         vol1 = volume(Zred1);
         %solution2 (axis-aligned):
         Zred2 = zonotope(interval(Z));
         Zred2 = repair(Zred2,Z);
         vol2 = volume(Zred2);

         if vol1<vol2
            P=polytope(Zred1, options);
         else
            P=polytope(Zred2, options);
         end
     end
% end

%repair zonotope if there is no length in one dimenion
function Zrep = repair(Z,Zorig)

%get length of each dimension
len = 2*rad(interval(Z));

%find zero lengths
index = find(len==0);

if ~isempty(index)
    %construct zonotope to be added
    origLen = 2*rad(interval(Zorig));
    origCenter = center(interval(Zorig));
    
    %get Zmatrix
    Zmat = Z.Z;
    
    for i=1:length(index)
        
        ind = index(i);
        %replace center
        Zmat(ind,1) = origCenter(ind);
        %replace generator value
        Zmat(ind,ind+1) = 0.5*origLen(ind);
    end

    %instantiate zonotopes
    Zrep = zonotope(Zmat);
else
    Zrep = Z;
end


%------------- END OF CODE --------------