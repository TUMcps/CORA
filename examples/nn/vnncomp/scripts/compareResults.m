% Specify the paths to the results files to compare.
resultsPath1 = 'results/250627-080921/cora';
resultsPath2 = '../zonotack-paper/neurips25evaluation/vnncomp2024_results/alpha_beta_crown';
% resultsPath2 = 'vnncomp2022_results/alpha_beta_crown';

% Find all benchmark result directories.
benchdirs = dir(resultsPath1);

for i=1:length(benchdirs)
    % Extract the benchmark name.
    benchname = benchdirs(i).name;

    % Read the first results file.
    try
        results1 = readtable(sprintf('%s/%s/results.csv', ...
            resultsPath1,benchname),'Delimiter',',');
        results1.Properties.VariableNames = ...
            {'name','onnxPath','vnnlibPath','prepTime','result','verifTime'};
    catch e
        fprintf('--- Cannot read results for "%s/%s"!\n', ...
            resultsPath1,benchname);
        continue;
    end
    
    if endsWith(results1{end,'onnxPath'},'test_nano.onnx')
        results1(end,:) = [];
    end
    
    % Read the second results file.
    try
        results2 = readtable(sprintf('%s/%s/results.csv', ...
            resultsPath2,benchname),'Delimiter',',');
        results2.Properties.VariableNames = ...
            {'name','onnxPath','vnnlibPath','prepTime','result','verifTime'};
    catch e
        fprintf('--- Cannot read results for "%s/%s"!\n', ...
            resultsPath2,benchname);
        continue;
    end
    
    if endsWith(results2{end,'onnxPath'},'test_nano.onnx')
        results2(end,:) = [];
    end

    fprintf('__________________________________________________________________\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf('BENCHMARK %s\n',benchname);
    fprintf('------------------------------------------------------------------\n');

    % Ensure that the results have the same number of entries.
    assert(height(results1) == height(results2));
    
    % Initialize results.
    numVerif = [0 0];
    numFals = [0 0];
    numUnknown = [0 0];

    for j=1:height(results1)
        % Extract the j-th row of the results.
        result1j = results1(j,:);
        result2j = results2(j,:);
        % Ensure we are comparing the same instance.
        assert(strcmp(result1j.name{1},result2j.name{1}));
        assert(strcmp(result1j.onnxPath{1},result2j.onnxPath{1}));
        assert(strcmp(result1j.vnnlibPath{1},result2j.vnnlibPath{1}));
        % Ensure there are no soundness issues.
        assert(~strcmp(result1j.result{1},'unsat') || ~strcmp(result2j.result{1},'sat'));
        assert(~strcmp(result1j.result{1},'sat') || ~strcmp(result2j.result{1},'unsat'));
    
        % Count results.
        numVerif = numVerif ...
            + [strcmp(result1j.result{1},'unsat') ...
                strcmp(result2j.result{1},'unsat')];
        numFals = numFals ...
            + [strcmp(result1j.result{1},'sat') ...
                strcmp(result2j.result{1},'sat')];
        numUnknown = numUnknown ...
            + [strcmp(result1j.result{1},'unknown') ...
                strcmp(result2j.result{1},'unknown')] ...
            + [strcmp(result1j.result{1},'timeout') ...
                strcmp(result2j.result{1},'timeout')];
    
        if strcmp(result1j.result{1},'unknown') ...
                & ~strcmp(result2j.result{1},'unknown') ...
                & ~strcmp(result2j.result{1},'timeout')
            fprintf('unknown instance (results 2: %s): %d\n', ...
                result2j.result{1},j);
        end
    end

    % Compute total number of instances.
    totalNum = numVerif + numFals + numUnknown;
    
    % Print summary for first results.
    table = CORAtableParameters(sprintf('RESULTS 1 (%s)',resultsPath1));
    table.printHeader();
    table.printContentRow('#Verified',sprintf('%d/%d [%.1f%%]', ...
        numVerif(1),totalNum(1),numVerif(1)/totalNum(1)*100));
    table.printContentRow('#Falsified',sprintf('%d/%d [%.1f%%]', ...
        numFals(1),totalNum(1),numFals(1)/totalNum(1)*100));
    table.printContentRow('#Unknown',sprintf('%d/%d [%.1f%%]', ...
        numUnknown(1),totalNum(1),numUnknown(1)/totalNum(1)*100));
    table.printMidBoundaryRow();
    table.printContentRow('Solved',sprintf('%d/%d [%.1f%%]', ...
        numVerif(1)+numFals(1),totalNum(1), ...
        (numVerif(1)+numFals(1))/totalNum(1)*100));
    table.printFooter();
    
    % Print summary for second results.
    table = CORAtableParameters(sprintf('RESULTS 2 (%s)',resultsPath2));
    table.printHeader();
    table.printContentRow('#Verified',sprintf('%d/%d [%.1f%%]', ...
        numVerif(2),totalNum(2),numVerif(2)/totalNum(2)*100));
    table.printContentRow('#Falsified',sprintf('%d/%d [%.1f%%]', ...
        numFals(2),totalNum(2),numFals(2)/totalNum(2)*100));
    table.printContentRow('#Unknown',sprintf('%d/%d [%.1f%%]', ...
        numUnknown(2),totalNum(2),numUnknown(2)/totalNum(2)*100));
    table.printMidBoundaryRow();
    table.printContentRow('Solved',sprintf('%d/%d [%.1f%%]', ...
        numVerif(2)+numFals(2),totalNum(2), ...
        (numVerif(2)+numFals(2))/totalNum(2)*100));
    table.printFooter();

end
