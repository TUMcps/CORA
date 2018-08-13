function [probTotal] = pyramid(pZ,mArray,P)
% pyramid - Encloses a probabilistic zonotope pZ by a pyramid with step
% sizes defined by an array of mSigma bounds and determines the probability
% of intersection with a polytope P
%
% Syntax:  
%    [pZ] = pyramid(pZ,mArray)
%
% Inputs:
%    pZ - probabilistic zonotope object
%    mArray - array of m-values for mSigma bounds
%    P - polytope object
%
% Outputs:
%    pZ - probabilistic zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 06-October-2007
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%get Sigma, dim
Sigma=sigma(pZ);
dim=length(Sigma);


%obtain array of max-values on mSigma bounds
for i=1:length(mArray)
    maxVal(i)=(2*pi)^(-0.5*dim)*det(Sigma)^(-0.5)*exp(-0.5*mArray(i)^2);
end
maxVal(end+1)=(2*pi)^(-0.5*dim)*det(Sigma)^(-0.5);

%obtain mSigma zonotopes
for i=1:length(mArray)
    msZ{i}=zonotope(pZ,mArray(i));
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
%     Znew=get(msZ{i},'Z');
%     Znew(3,1)=0.5*(maxVal(i+1)+maxVal(i));
%     Znew(3,end+1)=0.5*(maxVal(i+1)-maxVal(i));
%     Znew=zonotope(Znew);
%     V=vertices(Znew);
%     pP=polytope(get(V,'V')');
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
% IH=intervalhull([-10,10;-5,-3;-10,100]);
% P=polytope(IH);
% for i=1:length(mArray)
%     Znew=get(msZ{i},'Z');
%     Znew(3,1)=0.5*(maxVal(i+1)+maxVal(i));
%     Znew(3,end+1)=0.5*(maxVal(i+1)-maxVal(i));
%     Znew=zonotope(Znew);
%     V=vertices(Znew);
%     pP=polytope(get(V,'V')');
%     Pint=pP&P;
%     plot(Pint)
%     %plot3d(V);
%     hold on
% end

%------------- END OF CODE --------------