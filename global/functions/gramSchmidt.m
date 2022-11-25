function Q = gramSchmidt(B)
% gramSchmidt - construct an orthonormal basis Q from the vectors in B
%
% Syntax:  
%    Q = gramSchmidt(B)
%
% Inputs:
%    B - matrix with columns beeing the initial vectors for the basis
%
% Outputs:
%    Q - matrix whose columns span the orthonormal basis
%
% Example:
%    
%    B = [1 2;-1 1]
%    Q = gramSchmidt(B)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      23-January-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE -------------

    % construct suitable basis matrix if only a vector is provided
    if size(B,2) == 1
       B_ = eye(length(B));
       [~,ind] = max(abs(B'*B_));
       B_(:,ind) = [];
       B = [B,B_];
    end

    % gram schmidt method
    Q = zeros(size(B));
    R = zeros(size(B));

    for j = 1:size(B,2)
        v = B(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)'*B(:,j);
            v = v-R(i,j)*Q(:,i);
        end
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
 
end

%------------- END OF CODE --------------