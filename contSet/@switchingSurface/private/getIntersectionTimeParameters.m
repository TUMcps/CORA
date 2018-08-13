function [t_0, a_s, B_s, C_s] = getIntersectionTimeParameters(obj)
% getIntersectionTimeParameters - retrieve time model parameters from the y
% vector
%
% Syntax:  
%    [t_0, a_s, B_s, C_s] = getIntersectionTimeParameters(obj)
%
% Inputs:
%    obj - switchingSurface object
%
% Outputs:
%    t_0 - switching time constant
%    a_s - linear dependency on switching time
%    B_s - quadratic dependency on switching time
%    C_s - cubic dependency on switching time
%
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      22-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%init
t_0 = [];
a_s = [];
B_s = [];
C_s = [];

%write values to other variables
counter = 1;
t_0 = obj.y(counter);

%linear switching time parameter
for i = 1:obj.dim
    counter = counter + 1;
    a_s(i,1) = obj.y(counter);
end

if obj.timeOrder>1
    %quadratic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            counter = counter + 1;
            B_s(i,j) = obj.y(counter);
            %consider symmetric entries
            B_s(j,i) = obj.y(counter);
        end
    end  
end

if obj.timeOrder>2
    %cubic switching time parameter
    for i=1:obj.dim
        for j=i:obj.dim
            for k=j:obj.dim
                counter = counter + 1;
                C_s(i,j,k) = obj.y(counter);
                %consider symmetric entries
                C_s(i,k,j) = obj.y(counter);
                C_s(j,i,k) = obj.y(counter);
                C_s(j,k,i) = obj.y(counter);
                C_s(k,j,i) = obj.y(counter);
                C_s(k,i,j) = obj.y(counter);
            end
        end
    end
end
 

%------------- END OF CODE --------------