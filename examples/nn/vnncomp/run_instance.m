function [resStr,res] = run_instance(benchName,modelPath,vnnlibPath, ...
    resultsPath,timeout,verbose)
% run_instance - run the verification.
%
% Syntax:
%    [resStr,res] = run_instance(benchName,modelPath,vnnlibPath, ...
%       resultsPath,timeout,verbose)
%
% Inputs:
%    benchName - name of the benchmark
%    modelPath - path to the .onnx-file
%    vnnlibPath - path to the .vnnlib-file
%    resultsPath - path to the results directory
%    timeout - verification timeout
%
% Outputs:
%    resStr - result string
%    res - result code
%
% References:
%    [1] VNN-COMP'24
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       11-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    fprintf('run_instance(%s,%s,%s,%s,%d,%d)...\n',benchName,modelPath, ...
        vnnlibPath,resultsPath,timeout,verbose);

    % Initialize the result.
    res = struct('str','unknown','time',-1);
    % Count the number of verified branches.
    numVerified = 0;

    try
        % Measure verification time.
        totalTime = tic;

        fprintf('--- Loading MATLAB file...');
        % Create filename.
        [instanceFilename,modelName,~] = ...
            getInstanceFilename(benchName,modelPath,vnnlibPath);
        % Load stored network and specification.
        load(instanceFilename,'nn','options','permuteDims', ...
            'X0','specs');
        fprintf(' done\n');
        
        fprintf('--- Deleting MATLAB file...');
        % Delete file with stored networks and specification.
        delete(instanceFilename);
        fprintf(' done\n');

        % Obtain the model name.
        if permuteDims
            if strcmp(benchName,'collins_rul_cnn_2023') ...
                && ~strcmp(modelName,'NN_rul_full_window_40')
                inSize = nn.layers{1}.inputSize;
            else
                inSize = nn.layers{1}.inputSize([2 1 3]);
            end
            % inSize = nn.layers{1}.inputSize([2 1 3]);
        end

        fprintf('--- Running verification...');
        if verbose
            fprintf('\n\n');
        end

        % There can be multiple input sets. Concatenate the sets to a batch.
        x = [];
        r = [];
        for j=1:length(X0)
            % Extract the i-th input set.
            xi = 1/2*(X0{j}.sup + X0{j}.inf);
            ri = 1/2*(X0{j}.sup - X0{j}.inf);
            if permuteDims
                xi = reshape(permute(reshape(xi,inSize),[2 1 3]),[],1);
                ri = reshape(permute(reshape(ri,inSize),[2 1 3]),[],1);
            end
            % Append the i-th input set to the batch.
            x = [x xi];
            r = [r ri];
        end

        % Convert any safe set to a union of unsafe sets.
        specs = aux_convertSafeSet2UnionUnsafeSets(specs);

        if length(specs) > 1
            % We want to order the specification by criticallity. 
            % Therefore, we compute a quick adversarial attack for each 
            % specification.
            % Initialize an array for the criticallity scores.
            cv = zeros(length(specs),1);
            % Compute the sensitivity of the neural network.
            [S,~] = nn.calcSensitivity(x);
            % Compute the criticallity value for each specification.
            for i=1:length(specs)
                % Extract specification.
                [A,~,safeSet] = aux_spec2LinConstraint(specs(i));
                % Compute adversarial attacks.
                if safeSet
                    grad = pagemtimes(-A,S);
                else
                    grad = pagemtimes(A,S);
                end
                % Obtain number of constraints.
                [p,~] = size(A);
                % If there are multiple output constraints we try to falsify
                % each one individually.
                sgrad = reshape(permute(sign(grad),[2 3 1]),size(xi).*[1 p]);
        
                % Compute adversarial attacks based on the sensitivity.
                xi_ = repelem(xi,1,p) + repelem(ri,1,p).*sgrad;
                % We compute the outputs.
                yi_ = nn.evaluate(xi_);
                % We compute a criticallity score (< 0 means that the 
                % specification is violated).
                cv(i) = min(aux_computeCriticallityOfSpecs(specs(i),yi_));
            end
            % Sort the specification by their criticallity.
            [~,specIdx] = sort(cv,1,'ascend');
            specs = specs(specIdx);
        end

        % Keep track of unknown specifications.
        thereIsUnknown = false;

        % Handle multiple specs.
        for i=1:length(specs)
            % Extract specification.
            [A,b,safeSet] = aux_spec2LinConstraint(specs(i));
            
            while true
                try
                    % Reduce timeout in the case there are multiple 
                    % input sets.
                    remTimeout = timeout - toc(totalTime);
                    % Do verification.
                    [res,x_,y_] = nn.verify(x,r,A,b,safeSet, ...
                        options,remTimeout,verbose); % ,[2:3; 1:2]);
                    % Add the number of verified branches.
                    numVerified = numVerified + res.numVerified;
                    break;
                catch e
                    % Get the function name.
                    funcname = e.stack(1).name;
                    % Get the classname.
                    [dir,filename,~] = fileparts(e.stack(1).file);
                    if contains(dir,'@')
                        % The function is contained in a separate file.
                        [~,classname_] = fileparts(dir);
                        % Remove the '@'.
                        classname_(1) = [];
                        % Handle sub-functions.
                        if ~strcmp(filename,funcname)
                            % The error occurred in a sub-function.
                            funcname = [filename '/' funcname];
                        end
                        % Set the classname to the name of the parent 
                        % directory.
                        classname = classname_;
                    else
                        % The class name is the filename.
                        classname = filename;
                    end
                    % Get the line number.
                    linenr = e.stack(1).line;
                    % Print the error message. 
                    fprintf(newline);
                    fprintf( ...
                        'Unexpected Error! \n --- %s/%s [%d]: %s\n', ...
                        classname,funcname,linenr,e.message);
                    fprintf(newline);
                    % The verification failed.
                    res.str = 'UNKNOWN';
                    break;
                end
            end
            fprintf(' done (specification %d/%d)\n',i,length(specs));
            % Check if we could prove the unsafe set.
            switch res.str
                case 'VERIFIED'
                    % We have to verify the remaining specifications.
                    continue;
                case 'COUNTEREXAMPLE'
                    % We found a counterexample. We do not need to check the 
                    % remaining specifications.
                    thereIsUnknown = false;
                    break;
                otherwise
                    % We could not verify the specification. Continue
                    % trying to falsify other specifications.
                    thereIsUnknown = true;
            end
        end
    
        fprintf('Writing results...\n');
        fprintf('--- opening results file ...');
        % Open results file.
        fid = fopen(resultsPath,'w');
        fprintf(' done\n');

        if thereIsUnknown
            % There is an unknown specification.
            res.str = 'UNKNOWN';
        end
    
        fprintf('--- writing file ...');
        % Write results.
        switch res.str
            case 'VERIFIED'
                resStr = 'unsat';
                % Write content.
                fprintf(fid,['unsat' newline]);
                fclose(fid);
            case 'COUNTEREXAMPLE'
                resStr = 'sat';
                % Reorder input dimensions...
                if permuteDims
                  x_ = reshape(permute(reshape(x_,inSize([2 1 3])),[2 1 3]),[],1);
                end
                % Write content.
                fprintf(fid,['sat' newline '(']);
                % Write input values.
                for j=1:size(x_,1)
                    fprintf(fid,['(X_%d %f)' newline],j-1,x_(j));
                end
                % Write output values.
                for j=1:size(y_,1)
                    fprintf(fid,['(Y_%d %f)' newline],j-1,y_(j));
                end
                fprintf(fid,')');
                fclose(fid);
            otherwise
                resStr = 'unknown';
                % We cannot verify an input set; we dont have to check the other
                % input sets.
                fprintf(fid,['unknown' newline]);
                fclose(fid);
        end
        fprintf(' done\n');
    catch e
        fprintf(e.message);
        % There is some issue with the parsing; e.g. acasxu prop_6.vnnlib
        resStr = 'unknown';
        fprintf(' done\n');

        % Open results file.
        fid = fopen(resultsPath,'w');
        fprintf(fid,['unknown' newline]);
        fclose(fid);
    end

    if verbose
        % Print result.
        fprintf('%s -- %s: %s\n',modelPath,vnnlibPath,resStr);
        time = toc(totalTime);
        fprintf('--- Verification time: %.4f / %.4f [s]\n',time,timeout);
    end

    % Update the number of verified branches.
    res.numVerified = numVerified;
end


% Auxiliary functions -----------------------------------------------------

function cv = aux_computeCriticallityOfSpecs(specs,ys)
    % Check the specification and compute a value indicating how close we 
    % are to finding an adversarial example (< 0 mean the specification is 
    % violated).

    % Initialize the result.
    cv = inf([length(specs) 1]);

    for i=1:length(specs)
        % Extract the specification.
        [A,b,safeSet] = aux_spec2LinConstraint(specs(i));

        % Compute the logit difference.
        ld_ys = A*ys;
        % Compute the criticallity value per constraint.
        cvPerConstr = ld_ys - b;
        
        % Obtain the worst critical values, i.e., worst constraint, worst
        % input.
        if safeSet
            % safe iff all(A*y <= b) <--> unsafe iff any(A*y > b)
            % Thus, unsafe if any(-A*y < -b).
            cv(i) = min(-cvPerConstr,[],'all');
        else
            % unsafe iff all(A*y <= b) <--> safe iff any(A*y > b)
            % Thus, unsafe if all(A*y <= b).
            cv(i) = max(cvPerConstr,[],'all');
        end
    end
end

function unsafeUnionSpecs = aux_convertSafeSet2UnionUnsafeSets(specs)
    % Convert any safe set specifications to a union of unsafe sets.

    % Initialize the result.
    unsafeUnionSpecs = [];
    for i=1:length(specs)
        % Obtain the i-th specification.
        speci = specs(i);
        if strcmp(speci.type,'safeSet')
            % Obtain the i-th linear constraint.
            [A,b,~] = aux_spec2LinConstraint(speci);
            % We convert the safe set to unsafe sets, i.e., one for each
            % constraint.
            for j=1:size(A,1)
                % Invert the j-th constraint.
                Yj = polytope(-A(j,:),-b(j));
                % We add the inverted j-th constraint as an unsafe set.
                unsafeUnionSpecs = add(unsafeUnionSpecs,...
                    specification(Yj,'unsafeSet'));
            end
        else
            % We just add the unsafe set.
            unsafeUnionSpecs = add(unsafeUnionSpecs,speci);
        end
    end
end

function [A,b,safeSet] = aux_spec2LinConstraint(spec)
    % Extract a linear constraint from a specification.
    if isa(spec.set,'halfspace')
        A = spec.set.c';
        b = spec.set.d;
    else
        A = spec.set.A;
        b = spec.set.b;
    end
    % Obtain the type of specification.
    safeSet = strcmp(spec.type,'safeSet');
end

% ------------------------------ END OF CODE ------------------------------
