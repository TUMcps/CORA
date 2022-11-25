function res = testLongDuration_interval_mpower
% testLongDuration_interval_mpower - unit_test_function of power,
%    overloaded '^' operator for intervals
%
% Syntax:  
%    res = testLongDuration_interval_mpower
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
% See also: mtimes

% Author:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:      05-January-2016
% Last update:  08-August-2020 (MW, extend by random tests)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% Random tests ---------------------------------------------------------

numberRandTests = 10000;

% if exponent = 0, inf and sup are 1
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = a^0;
    if b.inf ~= 1 || b.sup ~= 1
        res = false;
        disp('test_mpower failed');
        return;
    end
end

% if exponent is even, inf is greater or equal to than 0
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    exponent = 2*randi([1,4],1,1); % 2,4,6,8
    b = a^exponent;
    if b.inf < 0
        res = false;
        disp('test_mpower failed');
        return;
    end
end

% if inf > 0, inf after operation also > 0
for i=1:numberRandTests
    randinf = randi(1000,1,1)*rand(1);
    a = interval(randinf,randinf*2);
    exponent = randi([1,5],1,1);
    b = a^exponent;
    if b.inf <= 0
        res = false;
        disp('test_mpower failed');
        return;
    end
end

% if exponent is non-integer, run through only if inf > 0
for i=1:numberRandTests
    randinf = randi(1000,1,1)*rand(1);
    a = interval(randinf,randinf*2);
    try
        b = a^rand(1);
        if b.inf <= 0
            res = false;
            disp('test_mpower failed');
            return;
        end
    catch
        res = false;
        disp('test_mpower failed');
        return;
    end
end

for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    try
        b = a^rand(1); % should throw an error
        % in case no error thrown, result should be NaN
        if ~all(isnan(b.inf))
            res = false;
            disp('test_mpower failed');
            return;
        end
    catch
        continue
    end
end


disp('test_mpower successful');
return;

%------------- END OF CODE --------------
