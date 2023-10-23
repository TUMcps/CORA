function res = test_interval_split
% test_interval_split - unit test function of split
%
% Syntax:
%    res = test_interval_split
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       29-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval
dim = 3;
lower = -rand(dim,1);
upper = rand(dim,1);
Int = interval(lower, upper);

% split interval along each dimension
for d=1:dim
    Int_split(d,:) = split(Int,d);
end

% true split
center = (upper + lower)/2;
lower_bottom = lower;
upper_bottom = [center(1); upper(2:3)];
lower_up = [center(1); lower(2:3)];
upper_up = upper;
Int_split_true{1,1} = interval(lower_bottom,upper_bottom);
Int_split_true{1,2} = interval(lower_up,upper_up);

upper_bottom = [upper(1); center(2); upper(3)];
lower_up = [lower(1); center(2); lower(3)];
Int_split_true{2,1} = interval(lower_bottom,upper_bottom);
Int_split_true{2,2} = interval(lower_up,upper_up);

upper_bottom = [upper(1:2); center(3)];
lower_up = [lower(1:2); center(3)];
Int_split_true{3,1} = interval(lower_bottom,upper_bottom);
Int_split_true{3,2} = interval(lower_up,upper_up);


% compare results
res = true;
for d=1:dim
    for s=1:2    % ...2 splits
        if Int_split{d,s} ~= Int_split_true{d,s}
            res = false;
            break
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
