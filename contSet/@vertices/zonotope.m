function [Z]=zonotope(obj)
% zonotopes - Computes a zonotope that encloses all vertices of the
% vertices object; the method is the one of principal component analysis
% developed by Olaf Stursberg in his 2003 HSCC paper
%
% Syntax:  
%    [Z]=zonotope(obj)
%
% Inputs:
%    obj - vertices object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    IH=intervalhull([1 2; -1 1]);
%    V=vertices(IH);
%    Z=zonotope(V);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author: Matthias Althoff
% Written: 09-May-2007
% Last update: 07-September-2007
% Last revision: ---

%------------- BEGIN CODE --------------

%if vertices object is nonempty
if ~isempty(obj.V)
    %compute the arithmetic mean of the vertices
    mean=sum(obj.V,2)/length(obj.V(1,:));

    %obtain sampling matrix
    translation=mean*ones(1,length(obj.V(1,:)));
    sampleMatrix=obj.V-translation;

    %compute the covariance matrix
    C=cov(sampleMatrix');

    %singular value decomposition
    [U,S,V] = svd(C);

    %auxiliary computations
    orientedMatrix=U'*sampleMatrix;
    m1=max(orientedMatrix,[],2);
    m2=min(orientedMatrix,[],2);

    %determine the center
    c=mean+U*(m1+m2)/2;


    %determine the generators
    for i=1:length(m1)
        G(:,i)=(m1(i)-m2(i))*0.5*U(:,i);
    end

    Z=zonotope([c,G]);
    
    %check...
    %Z2=parallelotope(obj,U);
else
    Z=[];
end

%------------- END OF CODE --------------