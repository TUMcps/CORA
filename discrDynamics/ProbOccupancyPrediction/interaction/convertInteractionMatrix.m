function [ThetaC]=convertInteractionMatrix(Theta)
% convertInteractionMatrix - converts the interaction matrix of a
% Markov-chain such that it can be used for an optimized update
%
% Syntax:  
%    [ThetaC]=convertInteractionMatrix(Theta)
%
% Inputs:
%    Theta - interaction matrix
%
% Outputs:
%    ThetaC - interaction matrix (optimized structure)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      12-October-2009
% Last update:  14-October-2009
%               10-August-2018
% Last revision:---

%------------- BEGIN CODE --------------

%obtain number of inputs and states
nrOfInputs_F=length(Theta(:,1));
nrOfInputs_L=length(Theta(1,:));
nrOfStates_F=length(Theta{1,1}(:,1));
nrOfStates_L=length(Theta{1,1}(1,:));

%allocate sparse Matrix and reserve space for nnz([Phi{:}]) entries
ThetaC=spalloc(nrOfInputs_F*nrOfStates_F,nrOfInputs_L*nrOfStates_L,nnz([Theta{:,:}]));

%reorganize structure
for iInput_F=1:nrOfInputs_F
    for iInput_L=1:nrOfInputs_L
        disp(['Input_F ',num2str(iInput_F),' Input_L ',num2str(iInput_L)]);

        %get nonzero entries
        [goalInd,initialInd]=find(Theta{iInput_F,iInput_L});

        %New row in the ThetaC matrix
        i=iInput_F+nrOfInputs_F*(goalInd-1);

        %New column in the ThetaC matrix
        j=iInput_L+nrOfInputs_L*(initialInd-1);

        %The entry that shall be copied to T(i,j) (as vector)
        s=Theta{iInput_F,iInput_L}(sub2ind(size(Theta{iInput_F,iInput_L}),goalInd,initialInd));

        ThetaC=ThetaC+sparse(i,j,s,nrOfInputs_F*nrOfStates_F,nrOfInputs_L*nrOfStates_L);
    end
end

%change ThetaC to full matrix
ThetaC=full(ThetaC);

%save results to file
[file,path] = uiputfile('*.mat','Save Interaction Matrix As');
cd(path);
save(file,'ThetaC');



%------------- END OF CODE --------------