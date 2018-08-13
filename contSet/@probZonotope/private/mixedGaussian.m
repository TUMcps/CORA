function fSum=mixedGaussian(x,omega,Sigma,mu)

dim=length(x);

for i=1:length(Sigma)
    f(i,:)=omega{i}/((2*pi)^(dim/2)*det(Sigma{i})^(1/2))*exp(-1/2*(x-mu{i})'*inv(Sigma{i})*(x-mu{i}));
end
fSum=sum(f);