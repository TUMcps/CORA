function [Y] = sample(E,N)
% sample - returns N samples that are contained within E
%
% Syntax:  
%    Y = sample(E,N) gives an array of samples which are contained in E
% Inputs:
%    E - ellipsoid object
%    N - #Samples
%
% Outputs:
%    Y - Array of sampled points
%
% Example: 
%    E = ellipsoid.generateRandom();
%    Y = sample(E,1000);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  16-October-2019 (complete revamp)
% Last revision:---

%------------- BEGIN CODE --------------
if E.isdegenerate
    [U,D] = eig(E.Q);
    d = diag(D);
    [d,ind] = sort(d,'descend');
    U = U(:,ind);
    n = size(E.Q,1);
    dred = d(1:end-(n-E.dim));
    %construct full-dimensional ellipsoid with smaller dimension
    %(x-c)'*inv(Q)*(x-c)<=1 -> y'*inv(L)*y<=1 with y=U'*(x-c)
    Et_red = ellipsoid(diag(dred));
    %sample reduced ellipsoid
    Yt_red = sample(Et_red,N);
    %add zeros (anything) to get back to original dimension
    Yt = [Yt_red;zeros(n-E.dim,size(Yt_red,2))];
    %back-transform to x-space
    Y = U*Yt + center(E);
else
    %computing box around E
    Z = zonotope(E);
    Y = zeros(E.dim,N);
    %compute the remaining N-nSamples samples
    nSamples = 0;
    p_acc = 1;
    sCounter = 0;
    %maximum amount of p_acc==0 loops
    maxZIt = 10;
    doInner = false;
    while nSamples<N
        %compute expected number of samples
        if p_acc == 0
            ns = N-nSamples;
        else
            ns = ceil((N-nSamples)*p_acc);
        end
        %compute samples for box
        X = sampleBox(Z,ns);
        %check which points are also in E
        bMask = containsPoint(E,X);
        bs = sum(bMask);
        p_acc = bs/ns;
        %this procedure only works if E is not too squished in any direction,
        %since then most of Z does not cover E
        %prevent infinite loop
        if p_acc==0
            if sCounter<maxZIt
                sCounter = sCounter + 1;
            else
                doInner = true;
                break;
            end
        end
        %save those points which are also in E
        Y(:,nSamples+1:nSamples+bs) = X(:,bMask);
        nSamples = bs + nSamples;
    end
    if doInner
        warning('Ellipsoid is too squished in one direction - samples are NOT uniformly distributed in E');
        %compute samples from the inner approximation
        Z_inner = zonotope(E,[],'u:box');
        Y = sampleBox(Z_inner,N);
    end
end
%------------- END OF CODE --------------