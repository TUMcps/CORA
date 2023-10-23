function res = test_spaceex2cora_parseForbidden()
% test_spaceex2cora_parseForbidden - example for parsing of forbidden
%    states/locations
%
% Syntax:
%    test_spaceex2cora_parseForbidden
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Maximilian Perschl
% Written:       22-September-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Since file reading etc. is tested in test_spaceex2cora_parseConfig,
% forbidden constellations are provided as strings and then given to the
% function to be tested

% linear constraint, no spec_map, equality
test_case(1) = "forbidden = ""x==10""";
solution{1}.specs = specification(polytope([],[],[1 0],100));
solution{1}.mapping = {};
state_names{1} = struct('name',{'x','v'});
component_names{1} = "comp1";
location_names{1} = {"loc1"};
error_expected(1) = false;
% linear constraint, no spec_map, inequality
test_case(2) = "forbidden = ""x<=10""";
solution{2}.specs = specification(polytope([1 0],10));
solution{2}.mapping = {};
state_names{2} = struct('name',{'x','v'});
component_names{2} = "comp1";
location_names{2} = {"loc1"};
error_expected(2) = false;
% linear constraint, no spec_map, inequality and equality
test_case(3) = "forbidden = ""x>=10&v==30""";
solution{3}.specs = specification(polytope([-1 0],-10,[0 1],30));
solution{3}.mapping = {};
state_names{3} = struct('name',{'x','v'});
component_names{3} = "comp1";
location_names{3} = {"loc1"};
error_expected(3) = false;
% linear constraint, with spec_map, inequality and equality
test_case(4) = "forbidden = ""x>=10&loc(comp1)=loc1|x==30&loc(comp1)=loc2""";
solution{4}.specs = specification(polytope([-1 0],-10));
solution{4}.specs(2) = specification(polytope([],[],[1 0],30));
solution{4}.mapping = {{2,3}};
state_names{4} = struct('name',{'x','v'});
component_names{4} = "comp1";
location_names{4} = {["loc1";"loc2"]};
error_expected(4) = false;
% nonlinear constraint, no spec_map, inequality
test_case(5) = "forbidden = ""x^2>=10""";
solution{5}.specs = specification(levelSet(str2sym("10-x^2"),sym(["x","v"]),"<="));
solution{5}.mapping = {};
state_names{5} = struct('name',{'x','v'});
component_names{5} = "comp1";
location_names{5} = {"loc1"};
error_expected(5) = true;


% test all cases
for i = 5:length(test_case)
    try
        [result_specs,result_mapping] = parseForbidden(test_case(i),...
            state_names{i},component_names{i},location_names{i});
        % compare result to expected result
        if ~( aux_compSpecs(result_specs(2:end),solution{i}.specs) ...
                && isequal(result_mapping,solution{i}.mapping) )
            res = false;
        end

    catch ME
        % check whether the test_case was supposed to throw an error
        if ~error_expected(i)
            res = false;
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

% Compare specification objects to interpret results
function equal = aux_compSpecs(spec1,spec2)
    equal = true;
    if length(spec1) ~= length(spec2)
        equal = false;
        return;
    end
    for i = 1:length(spec1)
        % Check set equality via two-way inclusion
        if ~(contains(spec2(i).set,spec1(i).set) && contains(spec1(i).set,spec2(i).set))
            equal = false;
            return;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
