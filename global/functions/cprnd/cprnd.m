function [X S] = cprnd(N,A,b,options)
%CPRND Draw from the uniform distribution over a convex polytope.
%   X = cprnd(N,A,b) Returns an N-by-P matrix of random vectors drawn
%   from the uniform distribution over the interior of the polytope
%   defined by Ax <= b. A is a M-by-P matrix of constraint equation
%   coefficients. b is a M-by-1 vector of constraint equation
%   constants.
%    
%   cprnd(N,A,b,options) allows various options to be specified.

%   'method'    Specifies the algorithm. One of the strings 'gibbs',
%               'hitandrun' (the default), and 'achr'. The default
%               algorithm 'hitandrun' is a vanilla hit-and-run sampler
%               [1]. 'gibbs' specifies a Gibbs sampler, 'achr' is the
%               Adaptive Centering Hit-and-Run algorithm of [2]. 
%
%   'x0'        Vector x0 is a starting point which should be interior 
%               to the polytope. If x0 is not supplied CPRND uses the
%               Chebyshev center as the initial point.
%    
%   'isotropic' Perform an adaptive istropic transformation of the 
%               polytope. The values 0, 1 and 2 respectively turn
%               off the transformation, construct it during a runup
%               period only, or continuously update the
%               tranformation throughout sample production. The
%               transformation makes sense only for the Gibbs and
%               Hit-and-run samplers (ACHR is invariant under
%               linear transformations).
%   
%   'discard'   Number of initial samples (post-runup) to discard.
%
%   'runup'     When the method is gibbs or hitandrun and the
%               isotropic transformation is used, this is the
%               number of initial iterations of the algorithm in
%               the untransformed space. That is a sample of size
%               runup is generated and its covariance used as the
%               basis of a transformation. 
%               When the method is achr runup is the number of
%               conventional hitandrun iterations. See [2].
%
%   'ortho'     Zero/one flag. If turned on direction vectors for
%               the hitandrun algorithm are generated in random
%               orthonormal blocks rather than one by
%               one. Experimental and of dubious utility.
%   
%   'quasi' Allows the user to specify a quasirandom number generator
%               (such as 'halton' or 'sobol').  Experimental and of
%               dubious utility.
%
%      

%   By default CPRND employs the hit-and-run sampler which may
%   exhibit marked serial correlation and very long convergence times.
%
% References
% [1] Kroese, D.P. and Taimre, T. and Botev, Z.I., "Handbook of Monte
% Carlo Methods" (2011), pp. 240-244.
% [2] Kaufman, David E. and Smith, Robert L., "Direction Choice for
% Accelerated Convergence in Hit-and-Run Sampling", Op. Res. 46,
% pp. 84-95.
%
% Copyright (c) 2011-2012 Tim J. Benham, School of Mathematics and Physics,
%                University of Queensland.

    function y = stdize(z)
        y = z/norm(z);
    end

    p = size(A,2);                         % dimension
    m = size(A,1);                         % num constraint ineqs
    x0 = [];
    runup = [];                              % runup necessary to method
    discard = [];                            % num initial pts discarded
    quasi = 0;
    method = 'achr';
    orthogonal = 0;
    isotropic = [];
    
    % gendir generates a random unit (direction) vector.
    gendir = @() stdize(randn(p,1));

    % Alternative function ogendir forces directions to come in
    % orthogonal bunches.
    Ucache = {};
    function u = ogendir()
        if length(Ucache) == 0
            u = stdize(randn(p,1));
            m = null(u');               % orthonormal basis for nullspace
            Ucache = mat2cell(m',ones(1,p-1));
        else
            u = Ucache{end}';
            Ucache(end) = [];
        end
    end
    
    % Check input arguments
    
    if m < p+1                             
        % obv a prob here
        error('cprnd:obvprob',['At least ',num2str(p+1),' inequalities ' ...
                            'required']);
    end

    if nargin == 4
        if isstruct(options)

            if isfield(options,'method')
                method = lower(options.method);
                switch method
                  case 'gibbs'
                  case 'hitandrun'
                  case 'achr'
                  otherwise
                    error('cprnd:badopt',...
                          ['The method option takes only the ' ...
                           'values "gibbs", "hitandrun", and "ACHR".']);
                end
            end
            
            if isfield(options,'isotropic')
                % Controls application of isotropic transformation,
                % which seems to help a lot.
                isotropic = options.isotropic;
            end

            if isfield(options,'discard')
                % num. of samples to discard
                discard = options.discard;
            end

            if isfield(options,'quasi')
                % Use quasi random numbers, which doesn't seem to
                % help much.
                quasi = options.quasi;
                if quasi && ~ischar(quasi), quasi='halton'; end
                if ~strcmp(quasi,'none')
                    qstr = qrandstream(quasi,p,'Skip',1);
                    gendir = @() stdize(norminv(qrand(qstr,1),0,1)');
                end
            end

            if isfield(options,'x0')
                % initial interior point
                x0 = options.x0;
            end

            if isfield(options,'runup')
                % number of points to generate before first output point
                runup = options.runup;
            end
            
            if isfield(options,'ortho')
                % Generate direction vectors in orthogonal
                % groups. Seems to help a little.
                orthogonal = options.ortho;
            end
            
        else
            x0 = options;                   % legacy support
        end
    end

    % Default and apply options
    
    if isempty(isotropic)
        if ~strcmp(method,'achr')
            isotropic = 2;
        else
            isotropic = 0;
        end
    end
    
    if orthogonal
        gendir = @() ogendir();
    end
    
    % Choose a starting point x0 if user did not provide one.
    if isempty(x0)
        x0 = chebycenter(A,b); % prob. if p==1?
    end
    
    % Default the runup to something arbitrary.
    if isempty(runup)
        if strcmp(method,'achr')
            runup = 10*(p+1);
        elseif isotropic > 0
            runup = 10*p*(p+1);
        else 
            runup = 0;
        end
    end

    % Default the discard to something arbitrary
    if isempty(discard)
        if strcmp(method,'achr')
            discard = 25*(p+1);
        else
            discard = runup;
        end
    end

    X = zeros(N+runup+discard,p);

    n = 0;                                  % num generated so far
    x = x0;

    % Initialize variables for keeping track of sample mean, covariance
    % and isotropic transform.
    M = zeros(p,1);                         % Incremental mean.
    S2 = zeros(p,p);                        % Incremental sum of
                                            % outer prodcts.
    S = eye(p); T = eye(p); W = A;
    
    while n < N+runup+discard
        y = x;
        
        % compute approximate stochastic transformation
        if isotropic>0
            if n == runup || (isotropic > 1 && n > runup)
                T = chol(S,'lower');
                W = A*T;
            end
            y = T\y;
        end
        
        switch method
            
          case 'gibbs'
            % choose p new components
            for i = 1:p
                % Find points where the line with the (p-1) components x_i
                % fixed intersects the bounding polytope.
                e = circshift(eye(p,1),i-1);
                z = W*e;
                c = (b - W*y)./z;
                tmin = max(c(z<0)); tmax = min(c(z>0));
                % choose a point on that line segment
                y = y + (tmin+(tmax-tmin)*rand)*e;
            end
            
          case 'hitandrun'
            % choose a direction
            u = gendir();
            % determine intersections of x + ut with the polytope
            z = W*u;
            c = (b - W*y)./z;
            tmin = max(c(z<0)); tmax = min(c(z>0));
            % choose a point on that line segment
            y = y + (tmin+(tmax-tmin)*rand)*u;
            
          case 'achr'
            % test whether in runup or not
            if n < runup
                % same as hitandrun
                u = gendir();
            else
                % choose a previous point at random
                v = X(randi(n),:)';
                % line sampling direction is from v to sample mean
                u = (v-M)/norm(v-M);
            end
            % proceed as in hit and run
            z = A*u;
            c = (b - A*y)./z;
            tmin = max(c(z<0)); tmax = min(c(z>0));
            % Choose a random point on that line segment
            y = y + (tmin+(tmax-tmin)*rand)*u;
        end
        
        if isotropic>0
            x = T*y;
        else
            x = y;
        end
        
        X(n+1,:)=x';
        n = n + 1;
        
        % Incremental mean and covariance updates
        delta0 = x - M;       % delta new point wrt old mean
        M = M + delta0/n;     % sample mean
        delta1 = x - M;       % delta new point wrt new mean
        if n > 1
            S2 = S2 + (n-1)/(n*n)*delta0*delta0'...
                 + delta1*delta1';
            S0 = S;
            S = S2/(n-1);           % sample covariance
        else 
            S = eye(p);
        end
        
    end
        
    X = X((discard+runup+1):(N+discard+runup),:); 

end
