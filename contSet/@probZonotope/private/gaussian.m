function f=gaussian(x,Sigma)

%get dimension
dim=length(Sigma);

%acceleration for 2-dim plots
if dim==2
    x1=x(1,:);
    x2=x(2,:);
    auxVar=1/det(Sigma)*(Sigma(2,2)*x1.^2-2*Sigma(1,2)*x1.*x2+Sigma(1,1)*x2.^2);
    f=1/((2*pi)^(dim/2)*det(Sigma)^(1/2))*exp(-1/2*auxVar);
%otherwise
else
    f=1/((2*pi)^(dim/2)*det(Sigma)^(1/2))*exp(-1/2*x'*inv(Sigma)*x);
end
