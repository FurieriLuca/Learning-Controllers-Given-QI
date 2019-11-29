function [U, cost_sample ] = sampling_ACC( cardinality, r, parameters,n,m,p,N,x0a,x0b,wa,wb,va,vb,C_b,P11,P12,M_b,R_b )
%Produces a cost sample

U = -1 + 2.*rand(cardinality,1);
while(norm(vec(U))<=0.01)
        U = -1 + 2.*rand(cardinality,1);
end
U=r*U/norm(vec(U));  %uniformly distributed random sample from the r-norm shell;


parameters_perturbed = parameters + U;


Kiterate = [parameters(1) 0 0 0 0 0;0 parameters(2) 0 0 0 0;0 0 parameters(1) 0 0 0;0 0 0 parameters(2) 0 0];


x0 = x0a+(x0b-x0a)*rand(n,1);
w_small = wa+(wb-wa)*rand(n*(N),1);
w=[x0;w_small];
v = va+(vb-va)*rand(p*(N+1),1);

y = C_b*inv(eye(size(P12,1))-P12*Kiterate*C_b)*P11*w+inv(eye(size(C_b,1))-C_b*P12*Kiterate)*v;
u = Kiterate*C_b*inv(eye(size(P12,1))-P12*Kiterate*C_b)*P11*w + Kiterate*inv(eye(size(C_b,1))-C_b*P12*Kiterate)*v;

cost_sample = y'*M_b*y + u'*R_b*u;


end

