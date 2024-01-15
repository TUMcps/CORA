function pass = test_spaceex2cora_parseInitial()
% test_spaceex2cora_parseInitial - example for parsing of initial sets
%
% Syntax:
%    test_spaceex2cora_parseInitial
%
% Inputs:
%    -
%
% Outputs:
%    pass - boolean

% Authors:       Maximilian Perschl
% Written:       22-September-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

pass = true;

% Since file reading etc. is tested in test_spaceex2cora_parseConfig,
% initial conditions are provided as strings and then given to the function
% to be tested

% basic case without inputs
test_case(1) = "initially = ""10<=x<=10.2 & v==0""";
solution{1}.states = interval([10,0],[10.2,0])';
solution{1}.inputs = interval.empty(1);
state_names{1} = struct('name',{'x','v'});
input_names{1} = [];
error_expected(1) = 0;
% basic case with inputs
test_case(2) = "initially = ""10<=x<=10.2 & v==0 & u == 10""";
solution{2}.states = interval([10,0],[10.2,0])';
solution{2}.inputs = interval(10,10);
state_names{2} = struct('name',{'x','v'});
input_names{2} = struct('name',{'u'});
error_expected(2) = 0;
% non fully specified case
test_case(3) = "initially = ""10<=x<=10.2""";
solution{3}.states = interval([10,0],[10.2,0])';
solution{3}.inputs = interval.empty(1);
state_names{3} = struct('name',{'x','v'});
input_names{3} = [];
error_expected(3) = 1;
% Multi line initially argument
test_case(4) = "initially = ""10<=x<=10.2 &" + newline + " v==0""";
solution{4}.states = interval([10,0],[10.2,0])';
solution{4}.inputs = interval.empty(1);
state_names{4} = struct('name',{'x','v'});
input_names{4} = [];
error_expected(4) = 0;
% Variable defined in 2 terms
test_case(5) = "initially = ""10<=x&x<=10.2""";
solution{5}.states = interval(10,10.2)';
solution{5}.inputs = interval.empty(1);
state_names{5} = struct('name',{'x'});
input_names{5} = [];
error_expected(5) = 0;


% test all cases
for i = 1:length(test_case)
    try
        [result_state,result_inputs] = parseInitial(test_case(i),0,state_names{i},input_names{i});
        % compare result to expected result
        if ~(isequal(result_state,solution{i}.states) && isequal(result_inputs,solution{i}.inputs))
            pass = false;
        end
    % if an error is thrown, check whether the test_case was supposed to
    % throw one
    catch
        if ~error_expected(i)
            pass = false;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
