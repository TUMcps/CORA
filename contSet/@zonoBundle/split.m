function Zsplit = split(zB,varargin)
% split - Splits a zonotope bundle into two zonotope bundles. This is done 
%    for one or every generator resulting in n possible splits where n is
%    the system dimension; it is also possible to use a splitting hyperplane
%
% Syntax:
%    Zsplit = split(zB,varargin)
%
% Inputs:
%    zB - zonoBundle object
%    N/hyperplane - splitting dimension in splitting hyperplane
%
% Outputs:
%    Zsplit - one or many zonotope bundle pairs
%
% Example: 
%    Z1 = zonotope([1;1;-1], [1 1 -1; -1 1 0; -2 0 1]);
%    Z2 = Z1 + [2;1;-1];
%    zB = zonoBundle({Z1,Z2}); 
%    Zsplit = split(zB,2);
%
%    figure; hold on;
%    plot(zB);
%    plot(Zsplit{1},[1,2],'LineStyle','--');
%    plot(Zsplit{2},[1,2],'LineStyle','--');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       03-February-2011
% Last update:   23-August-2013
%                25-January-2016
%                25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%split all dimensions
if nargin==1
    
    %obtain enclosing interval hull
    IH = interval(zB);
    %obtain limits
    leftLimit = infimum(IH);
    rightLimit = supremum(IH);
    
    for N = 1:length(leftLimit)
        %split one dimension of the interval hull
        Zsplit{N} = aux_splitOneDim(zB,leftLimit,rightLimit,N); 
    end
    
elseif nargin==2
    
    %split given dimension
    if isnumeric(varargin{1}) 
        N = varargin{1};
        %obtain enclosing interval hull
        IH = interval(zB);
        %obtain limits
        leftLimit = infimum(IH);
        rightLimit = supremum(IH);
        %split one dimension
        Zsplit = aux_splitOneDim(zB,leftLimit,rightLimit,N); 
        
    %split using a halfspace
    elseif isa(varargin{1},'halfspace')
        %obtain halfspace
        h = varargin{1};
        %obtain rotation matrix
        rotMat = aux_rotationMatrix(h);
        invRotMat = rotMat.';
        %obtain enclosing interval hull of rotated zonotope
        IH = interval(invRotMat*zB);
        intervals = get(IH,'intervals'); % dead: no asset property called 'intervals'
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
        Zsplit{1} = and_(zB,Znew{1},'exact');
        Zsplit{2} = and_(zB,Znew{2},'exact');
    end
    
elseif nargin==3
    
    if strcmp(varargin{2},'bundle')
        %split halfway in a direction using a zonotope bundle
        dir = varargin{1};
        Zsplit = aux_directionSplitBundle(zB,dir);
    end
    
end

end


% Auxiliary functions -----------------------------------------------------

function Zsplit = aux_splitOneDim(Zbundle,lb,ub,dim)

%split limits for a given dimension
leftLimitMod = lb;
leftLimitMod(dim) = 0.5*(lb(dim)+ub(dim));
rightLimitMod = ub;
rightLimitMod(dim) = 0.5*(lb(dim)+ub(dim));

%construct zonotopes which are the left and right boxes
Zleft = zonotope(interval(lb,rightLimitMod));
Zright = zonotope(interval(leftLimitMod,ub));

%generate splitted zonotope bundles
Zsplit{1} = and_(Zbundle,Zleft,'exact');
Zsplit{2} = and_(Zbundle,Zright,'exact');

%shrink zonotopes
%W{1} = eye(length(options.W));
% Zred = reduce(Zbundle.Z{1},'methC',1,options.filterLength);
% Zmat = Zred.Z;
% W{2} = Zmat(:,2:end); %<-- obtain this from order reduction!!

% Zsplit{1} = and_(Zbundle,shrink2(Zsplit{1},W),'exact');
% Zsplit{2} = and_(Zbundle,shrink2(Zsplit{2},W),'exact');

% Zsplit{1} = pinv(options.W)*shrink(options.W*Zsplit{1},options.filterLength);
% Zsplit{2} = pinv(options.W)*shrink(options.W*Zsplit{2},options.filterLength);

end

function Zsplit = aux_directionSplitBundle(Z,dir)

%obtain dimension
N = length(dir);

%obtain rotation matrix
newDir = [1; zeros(N-1,1)];
rotMat = aux_rotationMatrixDir(dir, newDir);

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
Zsplit{1} = zonoBundle(Z1);
Zsplit{2} = zonoBundle(Z2);

end

function rotMat = aux_rotationMatrix(h)

%get dimension
N = length(h.c);

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
    if N>2
        B(:,3:N) = null(B(:,1:2).'); 
    end
    
    %compute angle between uVec and n
    angle = acos(newDir.'*n);
    %rotation matrix
    R = eye(N);
    R(1,1) = cos(angle);
    R(1,2) = -sin(angle);
    R(2,1) = sin(angle);
    R(2,2) = cos(angle);
    %final rotation matrix
    rotMat = B*R*inv(B);
    
else
    if withinTol(h.c.'*newDir,1)
        rotMat = eye(N);
    else
        rotMat = -eye(N);
    end
end

end

function rotMat = aux_rotationMatrixDir(dir, newDir)

%get dimension
N = length(dir);

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
    if N>2
        B(:,3:N) = null(B(:,1:2).'); 
    end
    
    %compute angle between uVec and n
    angle = acos(newDir.'*n);
    %rotation matrix
    R = eye(N);
    R(1,1) = cos(angle);
    R(1,2) = -sin(angle);
    R(2,1) = sin(angle);
    R(2,2) = cos(angle);
    %final rotation matrix
    rotMat = B*R*inv(B);
    
else
     if withinTol(dir.'*newDir,1)
        rotMat = eye(N);
    else
        rotMat = -eye(N);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
