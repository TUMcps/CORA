function I = interval(SpS)
% interval - overapproximates a spectrahedral shadow by an interval
%
% Syntax:
%    I = interval(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    I - interval object
%
% Example: 
%    SpS = spectraShadow([eye(3) eye(3) eye(3)]);
%    I = interval(SpS);
%
%    figure; hold on;
%    plot(SpS,[1,2],'b');
%    plot(I,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       04-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isemptyobject(SpS) || representsa_(SpS,'emptySet',1e-9)
    I = interval.empty(dim(SpS));
    return
end

n = dim(SpS);

lower = zeros([n 1]);
upper = zeros([n 1]);

% Compute the support function in each axis direction
for i=1:n
    ei = zeros([n 1]);
    ei(i) = 1;
    
    lower(i) = SpS.supportFunc_(-ei,'upper');
    upper(i) = SpS.supportFunc_(ei,'upper');
end
I = interval(-lower, upper);

% ------------------------------ END OF CODE ------------------------------
