function cPZ = enclose(cPZ1,varargin)
% enclose - Generates a conPolyZono object that encloses a conPolyZono 
%           object and its linear transformation
%
% Syntax:  
%    cPZ = enclose(cPZ1,cPZ2)
%    cPZ = enclose(cPZ1,M,Splus)
%
% Inputs:
%    cPZ1 - first conPolyZono object
%    cPZ2 - second conPolyZono object, satisfying cPZ2 = (M*cPZ1) + cPZplus
%    M - matrix for the linear transformation
%    Splus - set added to the linear transformation
%
% Example: 
%    c = [0;0];
%    G = [1 0;0 1];
%    expMat = [1 0;0 1];
%    A = [1 -1];
%    b = 0;
%    expMat_ = [2 0;0 1];
%    cPZ1 = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    M = [1 2;-1 0]; b = [2;3];
%    cPZ2 = M*cPZ1 + b;
%   
%    cPZ = enclose(cPZ1,cPZ2);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor',[0.6 0.6 0.6],'Filled',true, ...
%         'EdgeColor','none','Splits',20);
%    plot(cPZ1,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%    plot(cPZ2,[1,2],'b','Filled',true,'EdgeColor','none','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/enclose, zonotope/enclose

% Author:       Niklas Kochdumper
% Written:      25-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin == 2
        cPZ2 = varargin{1};
    else
        M = varargin{1};
        Splus = varargin{2};
        if ~isa(Splus,'zonotope') && ~isa(Splus,'interval')
            [msg,id] = errWrongInput('Splus');
            error(id,msg);
        else
            cPZ2 = (M*cPZ1) + Splus;
        end
    end

    % check if exponent matrices are identical
    if ~all(size(cPZ1.id) == size(cPZ2.id)) ||  ...
       ~all(cPZ1.id == cPZ2.id) ||  ...
       ~all(size(cPZ1.expMat) == size(cPZ2.expMat)) || ...
       ~all(all(cPZ1.expMat == cPZ2.expMat)) || ...
       ~all(size(cPZ1.expMat_) == size(cPZ2.expMat_)) || ...
       ~all(all(cPZ1.expMat_ == cPZ2.expMat_)) || ...
       ~all(size(cPZ1.A) == size(cPZ2.A)) || ~all(all(cPZ1.A == cPZ2.A))

        error('Constraint polynomial zonotopes are not compatible!');
    end

    % compute set with the enclose method for polynomial zonotopes
    pZ1 = polyZonotope(cPZ1.c,cPZ1.G,cPZ1.Grest,cPZ1.expMat,cPZ1.id);
    pZ2 = polyZonotope(cPZ2.c,cPZ2.G,cPZ2.Grest,cPZ2.expMat,cPZ2.id);

    pZ = enclose(pZ1,pZ2);

    % construct resulting constrained polynomial zonotope
    expMat_ = [cPZ1.expMat_; zeros(1,size(cPZ1.expMat_,2))];
    
    cPZ = conPolyZono(pZ.c,pZ.G,pZ.expMat,cPZ1.A,cPZ1.b,expMat_, ...
                      pZ.Grest,pZ.id);
end

%------------- END OF CODE --------------