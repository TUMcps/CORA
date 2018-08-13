function [Zsplit] = split(Zbundle,varargin)
% split - Splits a zonotope bundle into two zonotope bundles. This is done 
% for one or every generator resulting in n possible splits where n is the 
% system dimension; it is also possible to use a splitting hyperplane
%
% Syntax:  
%    [Zsplit] = split(Zbundle,varargin)
%
% Inputs:
%    Z - zonotope bundle
%    dim/hyperplane - splitting dimension in splitting hyperplane
%    
%
% Outputs:
%    Zsplit - one or many zonotope bundle pairs
%
% Example: 
%    ---
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-February-2011
% Last update:  23-August-2013
%               25-January-2016
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%split all dimensions
if nargin==1
    %obtain enclosing interval hull
    IH = interval(Zbundle);
    %obtain limits
    leftLimit = infimum(IH);
    rightLimit = supremum(IH);
    
    for dim = 1:length(leftLimit)
        %split one dimension of the interval hull
        Zsplit{dim} = splitOneDim(Zbundle,leftLimit,rightLimit,dim); 
    end
elseif nargin==2
    %split given dimension
    if isnumeric(varargin{1}) 
        dim = varargin{1};
        %obtain enclosing interval hull
        IH = interval(Zbundle);
        %obtain limits
        leftLimit = infimum(IH);
        rightLimit = supremum(IH);
        %split one dimension
        Zsplit = splitOneDim(Zbundle,leftLimit,rightLimit,dim); 
    %split using a halfspace
    elseif isa(varargin{1},'halfspace')
        %obtain halfspace
        h = varargin{1};
        %obtain rotation matrix
        rotMat = rotationMatrix(h);
        invRotMat = rotMat.';
        %obtain enclosing interval hull of rotated zonotope
        IH = interval(invRotMat*Zbundle);
        intervals = get(IH,'intervals');
        %rotate halfspace
        h_rot = invRotMat*h;
        %new intervals
        newInterval{1} = intervals;
        newInterval{2} = intervals;
        newInterval{1}(1,:) = [intervals(1,1), h_rot.d];
        newInterval{2}(1,:) = [h_rot.d, intervals(1,2)]; 
        %new intervals
        intHull{1} = interval(newInterval{1});
        intHull{2} = interval(newInterval{2});
        %zonotope used for splitting
        Znew{1} = rotMat*zonotope(intHull{1});
        Znew{2} = rotMat*zonotope(intHull{2});
        %splitted sets
        Zsplit{1} = Zbundle & Znew{1};
        Zsplit{2} = Zbundle & Znew{2};
    end
elseif nargin==3
    if strcmp(varargin{2},'bundle')
        %split halfway in a direction using a zonotope bundle
        dir = varargin{1};
        Zsplit = directionSplitBundle(Zbundle,dir);
    end
end





function Zsplit = splitOneDim(Zbundle,leftLimit,rightLimit,dim)


%split limits for a given dimension
leftLimitMod = leftLimit;
leftLimitMod(dim) = 0.5*(leftLimit(dim)+rightLimit(dim));
rightLimitMod = rightLimit;
rightLimitMod(dim) = 0.5*(leftLimit(dim)+rightLimit(dim));

%construct zonotopes which are the left and right boxes
Zleft = zonotope(interval(leftLimit,rightLimitMod));
Zright = zonotope(interval(leftLimitMod,rightLimit));

%generate splitted zonotope bundles
Zsplit{1} = Zbundle & Zleft;
Zsplit{2} = Zbundle & Zright;   

%shrink zonotopes
%W{1} = eye(length(options.W));
% Zred = reduce(Zbundle.Z{1},'methC',1,options.filterLength);
% Zmat = get(Zred,'Z');
% W{2} = Zmat(:,2:end); %<-- obtain this from order reduction!!

% Zsplit{1} = Zbundle & shrink2(Zsplit{1},W);
% Zsplit{2} = Zbundle & shrink2(Zsplit{2},W);

% Zsplit{1} = pinv(options.W)*shrink(options.W*Zsplit{1},options.filterLength);
% Zsplit{2} = pinv(options.W)*shrink(options.W*Zsplit{2},options.filterLength);


function Zsplit = directionSplitBundle(Z,dir)

%obtain dimension
dim = length(dir);

%obtain rotation matrix
newDir = [1; zeros(dim-1,1)];
rotMat = rotationMatrixDir(dir, newDir);

%obtain enclosing interval
IH = interval(rotMat*Z);
intervals1 = get(IH,'intervals');
intervals2 = intervals1;

%split intervals
intervals1(1,2) = 0.5*(intervals1(1,1) + intervals1(1,2));
intervals2(1,1) = 0.5*(intervals2(1,1) + intervals2(1,2));
IH1 = interval(intervals1(:,1), intervals1(:,2));
IH2 = interval(intervals2(:,1), intervals2(:,2));

%zonotopes for zonotope bundle
Z1{1} = Z.Z{1};
Z1{2} = rotMat.'*zonotope(IH1);
Z2{1} = Z.Z{1};
Z2{2} = rotMat.'*zonotope(IH2);

%instantiate zonotope bundles
Zsplit{1} = zonotopeBundle(Z1);
Zsplit{2} = zonotopeBundle(Z2);


function rotMat = rotationMatrix(h)

%get dimension
dim = length(h.c);

if abs(h.c.'*newDir) ~= 1

    %normalize normal vectors
    n = h.c/norm(h.c);
    newDir = newDir/norm(newDir);
    %create mapping matrix
    B(:,1) = n;
    %find orthonormal basis for n, uVec
    indVec = newDir - (newDir.'*n)*n;
    B(:,2) = indVec/norm(indVec);
    %complete mapping matrix B
    if dim>2
        B(:,3:dim) = null(B(:,1:2).'); 
    end
    
    %compute angle between uVec and n
    angle = acos(newDir.'*n);
    %rotation matrix
    R = eye(dim);
    R(1,1) = cos(angle);
    R(1,2) = -sin(angle);
    R(2,1) = sin(angle);
    R(2,2) = cos(angle);
    %final rotation matrix
    rotMat = B*R*inv(B);
    
else
    if h.c.'*newDir == 1
        rotMat = eye(dim);
    else
        rotMat = -eye(dim);
    end
end

function rotMat = rotationMatrixDir(dir, newDir)

%get dimension
dim = length(dir);

if abs(dir.'*newDir) ~= 1

    %normalize normal vectors
    n = dir/norm(dir);
    newDir = newDir/norm(newDir);
    %create mapping matrix
    B(:,1) = n;
    %find orthonormal basis for n, uVec
    indVec = newDir - (newDir.'*n)*n;
    B(:,2) = indVec/norm(indVec);
    %complete mapping matrix B
    if dim>2
        B(:,3:dim) = null(B(:,1:2).'); 
    end
    
    %compute angle between uVec and n
    angle = acos(newDir.'*n);
    %rotation matrix
    R = eye(dim);
    R(1,1) = cos(angle);
    R(1,2) = -sin(angle);
    R(2,1) = sin(angle);
    R(2,2) = cos(angle);
    %final rotation matrix
    rotMat = B*R*inv(B);
    
else
    if dir.'*newDir == 1
        rotMat = eye(dim);
    else
        rotMat = -eye(dim);
    end
end


%------------- END OF CODE --------------