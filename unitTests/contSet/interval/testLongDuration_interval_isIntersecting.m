function res = testLongDuration_interval_isIntersecting
% testLongDuration_interval_isIntersecting - unit test function of
%    isIntersecting, note: only interval-to-interval tested!
%
% Syntax:  
%    res = testLongDuration_interval_isIntersecting
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Random cases
res_rand = true;
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi(50);

    % init random interval in [-1,1]
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % shift by center -> no intersection
    shift = ones(n,1);
    I_low = I - shift;
    I_high = I + shift;

    % check with correct solution
    if isIntersecting(I_low,I_high)
        res_rand = false; break;
    end

    % interval [-1,1]
    I = interval(-ones(n,1),ones(n,1));
    % shift by max. 1
    I_low = I - rand(1);
    I_high = I + rand(1);

    % check with correct solution
    if ~isIntersecting(I_low,I_high)
        res_rand = false; break;
    end
    
    % intervals meet at origin
    I = interval(-ones(n,1),ones(n,1));
    shift = ones(1,1);
    I_low = I - shift;
    I_high = I + shift;

    if ~isIntersecting(I_low,I_high)
        res_rand = false; break;
    end

end



% combine results
res = res_rand;

if res
    disp('test_isIntersecting successful');
else
    disp('test_isIntersecting failed');
end

%------------- END OF CODE --------------


