function probTotal = pyramid(probZ,mArray,P)
% pyramid - Encloses a probabilistic zonotope pZ by a pyramid with
%    step sizes defined by an array of mSigma bounds and determines
%    the probability of intersection with a polytope P
%
% Syntax:
%    probTotal = pyramid(pZ,mArray,P)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    mArray - array of m-values for mSigma bounds
%    P - polytope object
%
% Outputs:
%    probTotal - probabilistic zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       06-October-2007
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get Sigma, dimension
Sigma=sigma(probZ);
d = dim(probZ);

%obtain array of max-values on mSigma bounds
for i=1:length(mArray)
    maxVal(i)=(2*pi)^(-0.5*d)*det(Sigma)^(-0.5)*exp(-0.5*mArray(i)^2);
end
maxVal(end+1)=(2*pi)^(-0.5*d)*det(Sigma)^(-0.5);

%obtain mSigma zonotopes
for i=1:length(mArray)
    msZ{i}=zonotope(probZ,mArray(i));
end

%compute intersection probabilities
probTotal=0;
for i=1:length(mArray)
    %convert zonotope to polytope
    msP=polytope(msZ{i});
%     if i==1
%         plot(msZ{i});
%         hold on
%     end
    %intersect msP with P
    Pint=msP&P;
%     if chebyball(Pint)~=0
%         disp('check');
%     end
    %compute volume of intersection
    V=modVolume(Pint);
    %compute partial probability of intersection
    probPartial=(maxVal(i+1)-maxVal(i))*V;
    %add to total probability
    probTotal=probTotal+probPartial;
end

% plot(P)
% hold on
% %plot pyramid
% for i=1:length(mArray)
%     Znew=msZ{i}.Z;
%     Znew(3,1)=0.5*(maxVal(i+1)+maxVal(i));
%     Znew(3,end+1)=0.5*(maxVal(i+1)-maxVal(i));
%     Znew=zonotope(Znew);
%     V=vertices(Znew);
%     pP=polytope(V.V)');
%     plot(pP)
%     %plot3d(V);
%     hold on
% end
% 
% figure
% plot(P)
% hold on
% %plot remaining pyramid: only for the special case of the HSCC08 paper
% %example
% I=interval([-10,10;-5,-3;-10,100]);
% P=polytope(I);
% for i=1:length(mArray)
%     Znew=msZ{i}.Z;
%     Znew(3,1)=0.5*(maxVal(i+1)+maxVal(i));
%     Znew(3,end+1)=0.5*(maxVal(i+1)-maxVal(i));
%     Znew=zonotope(Znew);
%     V=vertices(Znew);
%     pP=polytope(V.V');
%     Pint=pP&P;
%     plot(Pint)
%     %plot3d(V);
%     hold on
% end

% ------------------------------ END OF CODE ------------------------------
