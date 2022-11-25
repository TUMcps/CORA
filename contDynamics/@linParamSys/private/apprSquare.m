function [Asquare] = apprSquare(A)
% square - computes the over-approximate square of an uncertain matrix 
%
% Syntax:  
%    [Asquare] = apprSquare(A)
%
% Inputs:
%    A - interval matrix defined by a set of matrices A_i
%
% Outputs:
%    Asquare - resulting interval matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      04-January-2009 
% Last update:  13-September-2016
% Last revision:---

%------------- BEGIN CODE --------------

%number of uncertain parameters
nrOfParam=length(A)-1;

%first method: simplify problem by abstarcting to interval matrices
[Aint] = intervalMatrix(A);

%compute square:
Asquare=exactSquare(Aint);
%AsquareOld=exactSquareOld(Aint);

%check with simpler computation:
%Asquare2=Aint*Aint;

%compute terms with common parameter:
%init int1, Asum1
int1=interval(0,1);
Asum1=0*A{1};

for i=1:nrOfParam
    Asum1=Asum1+int1*A{i+1}^2;
end

%compute terms with non-common parameters
%get number of possible combinations
comb=combinator(nrOfParam,2,'c');
nrOfComb=length(comb(:,1));

%init int2, Asum2
int2=interval(-1,1);
Asum2=0*A{1};

for i=1:nrOfComb
    %get indices
    ind1=comb(i,1)+1;
    ind2=comb(i,2)+1;
    
    %get V, absolute value of V
    V=A{ind1}*A{ind2}+A{ind2}*A{ind1};
    Vabs=abs(V);
    
    %add to Asum2
    Asum2=Asum2+Vabs;
end
Asum2=int2*Asum2;

%compute Asquare3
Asquare3=A{1}^2+A{1}*(Aint-A{1})+(Aint-A{1})*A{1}+Asum1+Asum2;

%COMMENT:
%Asquare3 seems to perform even worse for a single parameter!!

%compute intersection to use tightest over-approximation of both routines
Asquare=intersection(Asquare,Asquare3);

%------------- END OF CODE --------------