function [qZ] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of an object
% with a quadZonotope
%
% Syntax:  
%    [Z] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - quadZonotope object or numerical vector
%    summand2 - quadZonotope object or numerical vector
%
% Outputs:
%    Z - quadZonotpe after Minkowsi addition
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      04-September-2012
% Last update:  24-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%Find a quadZonotope object
%Is summand1 a quadZonotope?
if isa(summand1,'quadZonotope')
    %initialize resulting zonotope
    qZ = summand1;
    %initialize other summand
    summand = summand2;
%Is summand2 a quadZonotope?    
elseif isa(summand2,'quadZonotope')
    %initialize resulting zonotope
    qZ = summand2;
    %initialize other summand
    summand = summand1;  
end

%Is summand a zonotope?
if isa(summand,'quadZonotope')
    %auxiliary parameters
    p1 = length(qZ.G(1,:));
    p2 = length(summand.G(1,:));
    dim = length(qZ.c);
    
    %Calculate minkowski sum
    qZ.c = qZ.c + summand.c;
    qZ.G(:,(end+1):(end+length(summand.G))) = summand.G;
    qZ.Gsquare(:,(end+1):(end+length(summand.Gsquare))) = summand.Gsquare;
    qZ.Grest(:,(end+1):(end+length(summand.Grest))) = summand.Grest;
    
    %Gquad
    newGens = binom(p1+p2,2);
    Gquad = zeros(dim,newGens);
    counter = 0;
    counter1 = 0;
    counter2 = 0;

    for i = 1 : (p1+p2-1)
        for j = (i + 1) : (p1+p2)
            counter = counter + 1;
            if (i<=p1) && (j<=p1)
                counter1 = counter1 + 1;
                Gquad(:,counter) =  qZ.Gquad(:,counter1);
            elseif (i>p1) && (j>p1)
                counter2 = counter2 + 1;
                Gquad(:,counter) =  summand.Gquad(:,counter2);
            end
        end
    end
    qZ.Gquad = Gquad;

%Is summand a zonotope?
elseif isa(summand,'zonotope')
    %get zonotope matrix
    Zmat = get(summand,'Z');
    %Calculate minkowski sum
    qZ.c = qZ.c+Zmat(:,1);
    qZ.Grest(:,(end+1):(end+length(Zmat(1,2:end)))) = Zmat(:,2:end);
    
%is summand a vector?
elseif isnumeric(summand)
    %Calculate minkowski sum
    qZ.c = qZ.c+summand;
    
%something else?    
else
    qZ=[];
    disp('this operation is not implemented');
end

%------------- END OF CODE --------------