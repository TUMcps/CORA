function Zsplit = split(Z,varargin)
% split - splits a zonotope into two or more enclosing zonotopes
%
% Syntax:
%    Zsplit = split(Z,varargin)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zsplit - cell array of parallelpipeds represented as zonotopes
%
% Example: 
%    Z = zonotope(rand(2,4));
%    Zsplit = split(Z);
%
%    figure; hold on;
%    plot(Z);
%    plot(Zsplit{1}{1});
%    plot(Zsplit{1}{2});
%    plot(Zsplit{2}{1});
%    plot(Zsplit{2}{2});
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-January-2008
% Last update:   09-October-2008
%                26-June-2009
%                02-February-2011
%                27-August-2013
%                11-September-2013
%                25-January-2016
%                25-July-2016 (intervalhull replaced by interval)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize
if nargin == 1
    %split all dimensions
    Z=zonotope(interval(Z));
    for N=1:length(Z.c)
        %split one dimension
        Zsplit{N}=aux_splitOneDim(Z,N); 
    end
elseif nargin==2
    %no splitting halfspace is passed 
    if isnumeric(varargin{1})
        if isscalar(varargin{1})
            N=varargin{1};
            Z=zonotope(interval(Z));
            Zsplit=aux_splitOneDim(Z,N);
        else
            %split halfway in a direction
            dir = varargin{1};
            Zsplit = aux_directionSplit(Z,dir);
        end
    %split according to halfspace 
    else
        %obtain halfspace
        h = varargin{1};
        Zsplit = aux_halfspaceSplit(Z,h);
    end
elseif nargin==3
    if strcmp(varargin{2},'bundle')
        %split halfway in a direction using a zonotope bundle
        dir = varargin{1};
        Zsplit = aux_directionSplitBundle(Z,dir);
    elseif isscalar(varargin{1})
        N = varargin{1};
        maxOrder = varargin{2}; 
        %Zpp=pinv(options.W)*reduce(options.W*Z,'methD',maxOrder);
        Z=reduceGirard(Z,maxOrder); %<-- changed
        Zsplit=aux_splitOneDim(Z,N);
    else
        %compute split in perpendicular direction
        origDir = varargin{1};
        auxDir = varargin{2};
        
        %perpendicular direction
        projVal = auxDir.'*origDir/norm(origDir);
        perpDir = auxDir -  projVal*origDir/norm(origDir);
        
        %split halfway in perpendicular direction
        Zsplit = aux_directionSplit(Z,perpDir);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function Zsplit = aux_splitOneDim(Z,dim)

%center and generator matrix
c=Z.c;
G=Z.G;

%compute centers of splitted parallelpiped
c1=c-G(:,dim)/2;
c2=c+G(:,dim)/2;

%compute new set of generators
Gnew=G;
Gnew(:,dim)=Gnew(:,dim)/2;

%generate splitted parallelpipeds
Zsplit{1}=zonotope([c1,Gnew]);
Zsplit{2}=zonotope([c2,Gnew]);   

end

function Zsplit = aux_directionSplit(Z,dir)


%center and generator matrix
c = Z.c;
G = Z.G;

%aligned generator
alignedVal = 0;
dirUnitVec = dir/norm(dir);

for i = 1:length(G(1,:))
    %aligned part
    alignedPart = dirUnitVec.'*G(:,i);
    %enlarge aligned generator
    alignedVal = alignedVal + abs(alignedPart);
    %update generator
    G(:,i) = G(:,i) - alignedPart*dirUnitVec;
end

%new generator
newGen = alignedVal*dirUnitVec;

%beta value of aligned generator
beta = 0;

%check if intersection is possible
if abs(beta) < 1

    %new generators and centers
    g_1 = 0.5*(beta + 1)*newGen;
    g_2 = 0.5*(beta - 1)*newGen;
    c_1 = c + beta*newGen - g_1;
    c_2 = c + beta*newGen - g_2;


    %update zonotope
    Zsplit{1} = zonotope([c_1, g_1, G]);
    Zsplit{2} = zonotope([c_2, g_2, G]);
    
else
    Zsplit = Z; 
end

end

function Zsplit = aux_directionSplitBundle(Z,dir)

%obtain dimension
N = length(dir);

%obtain rotation matrix
newDir = [1; zeros(N-1,1)];
rotMat = aux_rotationMatrix(dir, newDir);

%obtain enclosing interval
IH = interval(rotMat*Z);
intervals1 = get(IH,'intervals');
intervals2 = intervals1;

%split intervals
intervals1(1,2) = 0.5*(intervals1(1,1) + intervals1(1,2));
intervals2(1,1) = 0.5*(intervals2(1,1) + intervals2(1,2));
IH1 = interval(intervals1);
IH2 = interval(intervals2);

%zonotopes for zonotope bundle
Z1{1} = Z;
Z1{2} = rotMat.'*zonotope(IH1);
Z2{1} = Z;
Z2{2} = rotMat.'*zonotope(IH2);

%instantiate zonotope bundles
Zsplit{1} = zonoBundle(Z1);
Zsplit{2} = zonoBundle(Z2);

end

function Zsplit = aux_halfspaceSplit(Z,h)

%halfspace values
dir = h.c;
d = h.d;

%center and generator matrix
c = Z.c;
G = Z.G;

%aligned generator
alignedVal = 0;
dirUnitVec = dir/norm(dir);

for i = 1:length(G(1,:))
    %aligned part
    alignedPart = dirUnitVec.'*G(:,i);
    %enlarge aligned generator
    alignedVal = alignedVal + abs(alignedPart);
    %update generator
    G(:,i) = G(:,i) - alignedPart*dirUnitVec;
end

%new generator
newGen = alignedVal*dirUnitVec;

%beta value of aligned generator
beta = (d - dir.'*c)/(dir.'*newGen);

%check if intersection is possible
if abs(beta) < 1

    %new generators and centers
    g_1 = 0.5*(beta + 1)*newGen;
    g_2 = 0.5*(beta - 1)*newGen;
    c_1 = c + beta*newGen - g_1;
    c_2 = c + beta*newGen - g_2;


    %update zonotope
    Zsplit{1} = zonotope([c_1, g_1, G]);
    Zsplit{2} = zonotope([c_2, g_2, G]);
    
else
    Zsplit = Z; 
end

end


function rotMat = aux_rotationMatrix(dir, newDir)

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
    if withinTol(dir.'*newDir,1)  %if dir.'*newDir == 1
        rotMat = eye(N);
    else
        rotMat = -eye(N);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
