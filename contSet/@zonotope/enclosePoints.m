function Z = enclosePoints(points,varargin)
% enclosePoints - enclose a point cloud with a zonotope
%
% Syntax:
%    Z = zonotope.enclosePoints(points)
%    Z = zonotope.enclosePoints(points,method)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%    method - method used for calculation
%               - 'maiga' (default)
%               - 'stursberg'
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    points = -1 + 2*rand(2,10);
%
%    Z1 = zonotope.enclosePoints(points);
%    Z2 = zonotope.enclosePoints(points,'maiga');
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'b');
%
% References:
%    [1] O. Stursberg et al. "Efficient representation and computation of 
%        reachable sets for hybrid systems", HSCC 2003
%    [2] M. Maiga et al. "A Comprehensive Method for Reachability Analysis
%        of Uncertain Nonlinear Hybrid Systems", TAC 2017
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       05-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    % parse input arguments
    method = setDefaultValues({'maiga'},varargin);

    % check input arguments
    inputArgsCheck({{points,'att','numeric','nonempty'};
                    {method,'str',{'maiga','stursberg'}}});

    % compute enclosing zonotope with the selected method
    if strcmp(method,'stursberg')
        Z = aux_enclosePointsStursberg(points);
    elseif strcmp(method,'maiga')
        Z = aux_enclosePointsMaiga(points);
    end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_enclosePointsStursberg(points)
% Computes an enclosing zonotope using the method from [1]

    % compute the arithmetic mean of the points
    mean = sum(points,2)/length(points(1,:));

    % obtain sampling matrix
    translation = mean*ones(1,length(points(1,:)));
    sampleMatrix = points-translation;

    % compute the covariance matrix
    C = cov(sampleMatrix');

    % singular value decomposition
    [U,~,~] = svd(C);

    % auxiliary computations
    orientedMatrix=U'*sampleMatrix;
    m1=max(orientedMatrix,[],2);
    m2=min(orientedMatrix,[],2);

    % determine the center
    c=mean+U*(m1+m2)/2;


    % determine the generators
    for i=1:length(m1)
        G(:,i)=(m1(i)-m2(i))*0.5*U(:,i);
    end

    Z = zonotope([c,G]);
end

function Zopt = aux_enclosePointsMaiga(points)
% Computes an enclosing zonotope using the method from [2]

    % initialization
    minVol = Inf;
    Zopt = [];
    N = 100;

    % loop over all sampling points
    for i = 1:N

        % compute current ratio
        r = i/N;

        % compute enclosing zonotope
        Z = aux_cloud2zonotope(points,r,size(points,1));

        % estimate the volume of the zonotope using the trace
        G = Z.G;
        tr = trace(G' * G);

        % update the minimum value
        if tr < minVol
            Zopt = Z;
            minVol = tr;
        end
    end
end

function Z = aux_cloud2zonotope(X,ratio,s1)
% implementation of the cloud2zonotope function in [2]

    n = size(X,1);
    c = zeros(n,1);
    iter = 0;
    R = [];
    w = zeros(size(X,1));
    r = zeros(size(X,1));
    while (iter < s1) 
        %|| condition(s2,w,rad))
        iter = iter + 1;
        I = interval.enclosePoints(X);
        mid = center(I);
        r = rad(I);
        X = X - mid;
        if iter == 1
           w = r; 
        end
        [U,~] = svd(X);
        u = U(:,1);
        g = ratio * abs(dot(u , r)) * u;
        X = aux_compress(X,g);
        c = c + mid;
        R = [R g];
    end
    I = interval.enclosePoints(X);
    mid = center(I);
    r = rad(I);
    c = c + mid;
    R = [R ,diag(r)];
    Z = zonotope([c,R]);
end

function X = aux_compress(X,g)
% implementation of the compress function in [2]
    
    u = g/norm(g);
    for i=1:size(X,2)
       d = dot(X(:,i) , u);
       d = min(norm(g),max(-norm(g),d));
       X(:,i) = X(:,i) - (d * u);
    end
end


% ------------------------------ END OF CODE ------------------------------
