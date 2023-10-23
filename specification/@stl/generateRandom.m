function obj = generateRandom(varargin)
% generateRandom - Generates a random STL-formula
%
% Syntax:
%    obj = stl.generateRandom()
%    obj = stl.generateRandom('Dimension',n)
%    obj = stl.generateRandom('Dimension',n,'NrOperators',nrOps)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension of the state space
%       <'NrOperators,nrOps> - number of logic opeartors ('&','next', etc.)
%       <'NestedOps',nrNest> - number of nested temporal operators
%       <'FinalTime',tFinal> - final time for the temporal logic formula 
%       <'TimeStep',dt> - time step size for the time specifications
%       <'Domain',dom> - domain in which the predicates should switch
%                       between true and false (class interval)
%
% Outputs:
%    obj - random STL-formula (class stl)
%
% Example: 
%    eq1 = stl.generateRandom()
%    eq2 = stl.generateRandom('Dimension',3)
%    eq3 = stl.generateRandom('Domain',interval([-1;-1],[1;1]))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       15-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % name-value pairs -> number of input arguments is a multiple of 2
    if mod(nargin,2) ~= 0
        throw(CORAerror('CORA:evenNumberInputArgs'));
    else
        % read input arguments
        NVpairs = varargin(1:end);
        % check list of name-value pairs
        checkNameValuePairs(NVpairs,{'Dimension','NrOperators','NestedOps', ...
                                        'FinalTime','TimeStep','Domain'});
        % dimension given?
        [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
        % number of opearators given?
        [NVpairs,nrOps] = readNameValuePair(NVpairs,'NrOperators');
        % number of nested opearators given?
        [NVpairs,nrNest] = readNameValuePair(NVpairs,'NestedOps');
        % final time given?
        [NVpairs,tFinal] = readNameValuePair(NVpairs,'FinalTime');
        % time step size given?
        [NVpairs,dt] = readNameValuePair(NVpairs,'TimeStep');
        % domain  given?
        [~,dom] = readNameValuePair(NVpairs,'Domain');
    end
    
    % default computation for dimension
    if isempty(n)
        if isempty(dom)
            nmax = 10;
            n = randi(nmax);
        else
            n = length(dom);
        end
    end
    
    % default number of operators
    if isempty(nrOps)
        opMax = 10;
        nrOps = randi(opMax);
    end

    % default final time
    if isempty(tFinal)
        tFinalMax = 10;
        tFinal = tFinalMax * rand();
    end

    % generate predicates
    list = generateRandomPredicates(n,nrOps,dom);

    % binary operations required to combine all predicates in the list
    nrBin = ceil(log2(nrOps));     

    % combine predicates using temporal operators
    cntUn = 0;
    unaryOps = {'next','globally','finally','~'};
    binaryOps = {'&','|','until','release'};
    timeNext = tFinal;

    for i = 1:nrOps

         % unary operation
        if length(list) == 1 || (rand() > 0.5 && cntUn < nrOps - nrBin)  

            ind = randi(length(list));
            
            if isempty(nrNest) || aux_nestedOps(list{ind}) < nrNest
                op = unaryOps{randi(length(unaryOps))};
            else
                op = unaryOps{randi(1)};
            end

            switch op
                case '~'
                    list{ind} = ~list{ind};
                case 'next'
                    time = min(timeNext,aux_randomTime(tFinal,dt));
                    timeNext = max(0,timeNext - time);
                    list{ind} = next(list{ind},time);
                case 'globally'
                    time = aux_randomTimeInterval(tFinal,dt);
                    timeNext = max(0,timeNext - supremum(time));
                    list{ind} = globally(list{ind},time);
                case 'finally'
                    time = aux_randomTimeInterval(tFinal,dt);
                    timeNext = max(0,timeNext - supremum(time));
                    list{ind} = finally(list{ind},time);
            end

            cntUn = cntUn + 1;

        % binary operation
        else                                        

            ind1 = randi(length(list));
            ind2 = randi(length(list));

            if isempty(nrNest) || ...
                   aux_nestedOps(list{ind1}) + aux_nestedOps(list{ind2}) < nrNest
                op = binaryOps{randi(length(binaryOps))};
            else
                op = binaryOps{randi(2)};
            end

            switch op
                case '&'
                    list{ind1} = list{ind1} & list{ind2};
                    list{ind2} = [];
                case '|'
                    list{ind1} = list{ind1} | list{ind2};
                    list{ind2} = [];
                case 'until'
                    time = aux_randomTimeInterval(tFinal,dt);
                    timeNext = max(0,timeNext - supremum(time));
                    list{ind1} = until(list{ind1},list{ind2},time);
                    list{ind2} = [];
                case 'release'
                    time = aux_randomTimeInterval(tFinal,dt);
                    timeNext = max(0,timeNext - supremum(time));
                    list{ind1} = until(list{ind1},list{ind2},time);
                    list{ind2} = [];
            end

            list = list(~cellfun('isempty',list));
        end       
    end

    obj = list{1};
end


% Auxiliary functions -----------------------------------------------------

function res = aux_randomTime(tFinal,dt)
% generate random time

    if isempty(dt)
        res = tFinal * rand();
    else
        res = max(0,dt * randi(floor(tFinal/dt-1)));
    end
end

function res = aux_randomTimeInterval(tFinal,dt)
% generate random time interval

    time1 = aux_randomTime(tFinal,dt);
    time2 = aux_randomTime(tFinal,dt);

    while time1 == time2
        time2 = aux_randomTime(tFinal,dt);
    end

    res = interval(min([time1,time2]),max([time1,time2]));
end

function res = aux_nestedOps(obj)
% determine the number of nested temporal operators
    
    if ~obj.temporal

        res = 0;

    elseif strcmp(obj.type,'finally') || ...
            strcmp(obj.type,'globally') || strcmp(obj.type,'next')

        res = aux_nestedOps(obj.lhs) + 1;

    elseif strcmp(obj.type,'until') || strcmp(obj.type,'release')

        res = aux_nestedOps(obj.lhs) + aux_nestedOps(obj.rhs) + 1;

    elseif strcmp(obj.type,'~')

        res = aux_nestedOps(obj.lhs);

    elseif strcmp(obj.type,'&') || strcmp(obj.type,'|')

        res = aux_nestedOps(obj.lhs) + aux_nestedOps(obj.rhs);

    end
end

% ------------------------------ END OF CODE ------------------------------
