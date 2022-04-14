function I = minkDiff(I1,I2,varargin)
% minkDiff - compute the Minkowski difference of two intervals:
%            I1 - I2 = I <-> I + I2 \subseteq I1
%
% Syntax:  
%    I = minus(I1,I2)
%    I = minus(I1,I2,type)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'approx')
%
% Outputs:
%    I - interval object after Minkowski difference
%
% Example: 
%    I1 = interval([-2;-1],[3;3]);
%    I2 = interval([-1;-1],[1;1]);
%
%    I = minkDiff(I1,I2);
%
%    figure; hold on;
%    plot(I1);
%    plot(I2,[1,2],'r');
%    plot(I,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minus

% Author:       Niklas Kochdumper
% Written:      10-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % different algorithms for different set representations
    if isnumeric(I2)
        
       I = I1 + (-I2);
       
    elseif isa(I2,'interval')
        try
           I = interval(infimum(I1)-infimum(I2),supremum(I1)-supremum(I2)); 
        catch
           I = []; 
        end
        
    else
        
        % parse input arguments
        type = 'approx';
        if nargin > 2 && ~isempty(varargin{1})
            type = varargin{1};
        end
        
        % check input arguments
        if strcmp(type,'exact')
           error(errNoExactAlg(I1,I2));
        end
        
        % compute inner-approximation of the Minkowski difference
        infi = infimum(I1); sup = supremum(I1); n = dim(I1);
        
        for i = 1:n
           temp = zeros(n,1);
           temp(i) = 1;
           sup(i) = sup(i) - supportFunc(I2,temp,'upper');
           infi(i) = infi(i) + supportFunc(I2,-temp,'upper');
        end
        
        % construct resulting interval
        try
            I = interval(infi,sup);
        catch
            I = []; 
        end
    end
end

%------------- END OF CODE --------------