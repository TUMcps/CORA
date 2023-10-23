function P_out = generateRandom(varargin)
% generateRandom - Generates a random polytope
%    restrictions:
%    - number of constraints >= 1    =>  can be unbounded
%    - number of constraints >= 2    =>  + can be unbounded and degenerate
%    - number of constraints >= n+1  =>  + can be bounded
%    - number of constraints >= n+2  =>  + can be bounded and degenerate
%
% Syntax:
%    P_out = polytope.generateRandom()
%    P_out = polytope.generateRandom('Dimension',n)
%    P_out = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon)
%    P_out = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon,...
%           'IsDegenerate',isDeg)
%    P_out = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon,...
%           'IsDegenerate',isDeg,'IsBounded',isBounded)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NrConstraints',nrCon> - number of constraints
%       <'IsDegenerate',isDeg> - degeneracy (true/false)
%       <'IsBounded',isBounded> - boundedness (true/false)
%
% Outputs:
%    P_out - random polytope
%
% Example: 
%    P = polytope.generateRandom('Dimension',2);
%    plot(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       05-May-2020
% Last update:   10-November-2022 (MW, adapt to new name-value pair syntax)
%                30-November-2022 (MW, new algorithm)
%                12-December-2022 (MW, less randomness for more stability)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','NrConstraints','IsDegenerate','IsBounded'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension',...
        @(x) mod(x,1) == 0 && x > 0);
    % number of constraints given?
    [NVpairs,nrCon] = readNameValuePair(NVpairs,'NrConstraints',...
        @(x) mod(x,1) == 0 && x > 0);
    % degeneracy given?
    [NVpairs,isDeg] = readNameValuePair(NVpairs,'IsDegenerate','islogical');
    % boundedness given?
    [NVpairs,isBnd] = readNameValuePair(NVpairs,'IsBounded','islogical');
end

% quick exits
if ~isempty(nrCon)
    if nrCon == 1
        % only one constraint -> cannot be bounded or degenerate
        if ~isempty(isBnd) && isBnd
            throw(CORAerror('CORA:wrongValue','name-value pair isBounded',...
                'cannot be true if only one constraint given.'));
        elseif ~isempty(isDeg) && isDeg
            throw(CORAerror('CORA:wrongValue','name-value pair isDegenerate',...
                'cannot be true if only one constraint given.'));
        end
    end
    if ~isempty(n)
        % dimension given -> cannot be bounded unless nrCon >= n+1
        if ~isempty(isBnd) && isBnd && nrCon < n+1
            throw(CORAerror('CORA:wrongValue','name-value pair isBounded',...
                'cannot be true if number of constraints < dimension + 1'));
        end
    end
end
if ~isempty(n) && n == 1 && ~isempty(isBnd) && ~isBnd && ~isempty(isDeg) && isDeg
    % one-dimensional polytopes cannot be degenerate and unbounded
    throw(CORAerror('CORA:wrongValue','name-value pair isDegenerate/isBounded',...
        'one-dimensional polytopes cannot be unbounded and degenerate at the same time'));
end

% default computation for dimension
if isempty(n)
    if ~isempty(nrCon) && (isempty(isBnd) || isBnd)
        % number of constraints given and boundedness either true by
        % default or because set so by user)
        if ~isempty(isDeg) && isDeg
            % degenerate: n <= number of constraints - 1
            n = randi([1,nrCon-2]);
        else
            % non-degenerate (default): n <= number of constraints - 2
            n = randi([1,nrCon-1]);
        end
    else
        % either: number of constraints not given
        % or:     number of constraints given, but unbounded
        % -> dimension can be set to any random value
        maxdim = 10;
        if ~isempty(isBnd) && ~isBnd && ~isempty(isDeg) && isDeg
            % a set can only be unbounded and degenerate for n > 1
            n = 2;
%             n = randi([2,maxdim]);
        else
            n = randi(maxdim);
        end
    end
else
    % dimension is given
    if ~isempty(nrCon)
        % number of constraints given: decide boundedness
        if ~isempty(isBnd) && isBnd
            % user wants the polytope to be bounded
            if ~isempty(isDeg) && isDeg
                % user wants the polytope to be degenerate
                if nrCon < n+2
                    throw(CORAerror('CORA:wrongValue','name-value pair isDegenerate',...
                        'has to be false if the polytope should be bounded and the number of constraints < dimension + 2.'));
                end
            else
                % only bounded, degeneracy does not matter
                if nrCon < n+1
                    throw(CORAerror('CORA:wrongValue','name-value pair isBounded',...
                        'has to be false if number of constraints < dimension + 1.'));
                end
            end
            % ...isBounded = true;
        else
            % either boundedness not defined by user or set to false
            if isempty(isBnd)
                % set to true/false depending on nrCon < n+1?
                if nrCon < n+1
                    isBnd = false;
                else
                    isBnd = true;
                end
            elseif nrCon < n+1
                % ensure that nrCon >= n+1
                throw(CORAerror('CORA:wrongValue','name-value pair isBounded',...
                        'has to be false if number of constraints < dimension + 1.'));
            end
        end
    end
end

% default degeneracy
if isempty(isDeg)
    isDeg = false;
end

% default boundedness
if isempty(isBnd)
    isBnd = true;
end

% default computation for number of constraints
if isempty(nrCon)
    if n == 1 && isDeg
        nrCon = 2;
    else
        nrCon = 2*n + randi([floor(n/2),n]);
    end
else
    % number of constraints are given -> ensure that other provided values
    % are admissible
    if ~isempty(isBnd) && isBnd && nrCon < n+1
        throw(CORAerror('CORA:wrongValue','name-value pair IsBounded',...
            'cannot be true if number of constraints < n+1.'));
    end
end


if ~isBnd
    % idea for unbounded case: sample directions

    % degenerate case: one pair of constaints needs to be of the form
    %    a*x <= b   a*x >= b

    if n == 1 && ~isDeg
        % 1D algorithm: can only be bounded at one side (towards +/-Inf)
        s = sign(randn); % randn will practically never be 0...
        A = s*ones(nrCon,n);
        b = s*rand(nrCon,1);

        % note that all but one constraint will be redundant
        P_out = polytope(A,b);
        
        % set properties
        P_out.bounded.val = isBnd;
        P_out.fullDim.val = ~isDeg;
        P_out.emptySet.val = false;
        return

    elseif n == 2 && ~isDeg
        % 2D: sample directions from one half of the plane, offset > 0
        A = zeros(nrCon,n);
        b = zeros(nrCon,1);
        for i=1:nrCon
            A(i,:) = [randn rand];
            b(i) = rand;
        end

        % map by random matrix
        M = randn(n);

        % init polytope
        P_out = polytope(A*M^(-1),b);

        % set properties
        P_out.bounded.val = isBnd;
        P_out.fullDim.val = ~isDeg;
        P_out.emptySet.val = false;
        return

    elseif n == 2 && isDeg
        % 2D: unbounded and degenerate -> line (only allows for constraints
        % along the same parallel vector, fill up with redundancies until
        % number of constraints reached)

        % random vector
        vec = randn(n,1);
        % normalize
        vec = vec / vecnorm(vec);

        % random offset
        offset = randn(1);

        % init constraint matrix and offset
        A = [vec'; -vec'; zeros(nrCon-2,n)];
        b = [offset; -offset; zeros(nrCon-2,1)];
        
        % add redundant halfspace
        for i=1:nrCon-2
            s = sign(randn);
            A(2+i,:) = s*vec;
            b(2+i) = s*offset + rand;
        end

        % map by random matrix
        M = randn(n);

        % init polytope
        P_out = polytope(A*M^(-1),b);

        % set properties
        P_out.bounded.val = isBnd;
        P_out.fullDim.val = ~isDeg;
        P_out.emptySet.val = false;
        return

    else
        % all higher-dimensional cases: sample directions, but keep one
        % entry always 0 -> unboundness in one dimension ensured
        
        % init constraint matrix and offset in dimension n-1
        tempdirs = randn(nrCon,n-1);
        A = (tempdirs' ./ vecnorm(tempdirs'))';
        b = ones(nrCon,1);

        % randomly selected dimension
        randDim = randi(n);

        % set all entries in that dimension to 0
        A = [A(:,1:randDim-1), zeros(nrCon,1), A(:,randDim:end)];

        % degenerate case
        if isDeg
            % adapt last one using second-to-last one
            A(end-1,:) = -A(end,:);
            b(end-1) = -b(end);

            % permutate matrix, offset
            order = randperm(nrCon);
            A = A(order,:);
            b = b(order);
        end

        % map by random matrix
        M = randn(n);

        % init polytope
        P_out = polytope(A*M^(-1),b);

        % set properties
        P_out.bounded.val = isBnd;
        P_out.fullDim.val = ~isDeg;
        P_out.emptySet.val = false;
    end
    
else
    % bounded case: generate simplex and add constraints
    % it is ensured that nrCon > n+1 (non-degenerate)
    %                    nrCon > n+2 (degenerate)

    if n == 1 && isDeg
        % has to be a point; add redundant constraints until number of
        % constraints
        point = randn;
        A = [1; -1; ones(nrCon-2,n)];
        b = [point; -point; rand(nrCon-2,1)];

        % init polytope
        P_out = polytope(A,b);

        % set properties
        P_out.bounded.val = isBnd;
        P_out.fullDim.val = ~isDeg;
        P_out.emptySet.val = false;
        return
    end

    % step 1: simplex
    nrConSimplex = n+1;
    A = [eye(n); -ones(1,n)/sqrt(n)];
    b = ones(nrConSimplex,1);
    
    % additional number of constraints
    addCon = nrCon - nrConSimplex;
    A = [A; zeros(addCon,n)];
    b = [b; zeros(addCon,1)];
    
    % step 2: sample random directions of length 1, then the resulting
    %         halfspace is always non-redundant
    for i=1:(nrCon-nrConSimplex)
        % random direction
        tempdir = randn(1,n);
        A(nrConSimplex+i,:) = tempdir / vecnorm(tempdir);
        % set offset to 1 to ensure non-redundancy
        b(nrConSimplex+i) = 1;
    end

    % degenerate case
    if isDeg
        % permutate matrix, offset
        order = randperm(nrCon);
        A = A(order,:);
        b = b(order);

        % adapt last one using second-to-last one
        A(end-1,:) = -A(end,:);
        b(end-1) = -b(end);
    end

    % step 3: rotate
    [M,~,~] = svd(randn(n)); % invertible
    
    % instantiate polytope
    P_out = polytope(A*M^(-1),b);
    
    % step 4: translate polytope by some small offset
    z = randn(n,1);
    P_out = P_out + z;

    % set properties
    P_out.bounded.val = isBnd;
    P_out.fullDim.val = ~isDeg;
    P_out.emptySet.val = false;

end


% ------------------------------ END OF CODE ------------------------------
