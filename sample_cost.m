%%%sample a random trajectory

U = -1 + 2.*rand(cardinality,1);
while(norm(vec(U))<=0.01)
        U = -1 + 2.*rand(cardinality,1);
end
U=r*U/norm(vec(U));  %uniformly distributed random sample from the r-norm shell;


parameters_perturbed = parameters + U;


rebuild_K;


x0 = x0a+(x0b-x0a)*rand(n,1);
w_small = wa+(wb-wa)*rand(n*(N),1);
w=[x0;w_small];
v = va+(vb-va)*rand(p*(N),1);

x = inv(eye(size(P12,1))-P12*Kiterate*C_b)*(P11*w+P12*Kiterate*v);
u = Kiterate*C_b*inv(eye(size(P12,1))-P12*Kiterate*C_b)*P11*w + Kiterate*inv(eye(size(C_b,1))-C_b*P12*Kiterate)*v;

cost_sample = x'*M_b*x + u'*R_b*u;


