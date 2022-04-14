function C = plus(summand1,summand2)
% plus - Overloaded '+' operator for the over-approximative Minkowski 
% addition of two capsules or the exact translation of a capsule by a 
% vector
%
% Syntax:  
%    C = plus(summand1,summand2)
%
% Inputs:
%    summand1 - capsule or numerical vector
%    summand2 - capsule or numerical vector
%
% Outputs:
%    C - capsule after Minkowsi addition
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    summand1 = C;
%    summand2 = [2; 2; 1];
%    C1 = C + summand1;
%    C2 = C + summand2;
%    plot(C);
%    hold on
%    plot(C1);
%    plot(C2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      04-March-2019 
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

%Find a capsule object
%Is summand1 a capsule?
if isa(summand1,'capsule')
    %initialize resulting capsule
    C = summand1;
    %initialize other summand
    summand = summand2;
%Is summand2 a capsule?    
elseif isa(summand2,'capsule')
    %initialize resulting capsule
    C = summand2;
    %initialize other summand
    summand = summand1;  
end

%Is summand a capsule?
if isa(summand,'capsule')
    % add centers
    C.c = C.c + summand.c;
    % only replace generator if the one of the summand is longer
    length_C = norm(C.g);
    length_summand = norm(summand.g);
    if length_C < length_summand
        C.g = summand.g;
        radiusOfGenerator = length_C;
    else
        radiusOfGenerator = length_summand; 
    end
    % obtain new radius
    C.r = C.r + summand.r + radiusOfGenerator;
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate Minkowski sum
    C.c = C.c + summand;
    
elseif isa(summand,'conPolyZono')
    
    C = summand + C;
    
%something else?    
else
    % throw error for given arguments
    error(noops(summand1,summand2));
end

%------------- END OF CODE --------------