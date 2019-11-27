
%%%THIS IS A STANDALONE FILE that verifies gradient dominance beyond QI

clear all
clc;
A = [3 -2;-2 -1];

n=size(A,1);
B = [-5 -1;2 5];

C = eye(n);

m=size(B,2);
p=size(C,1);

N=2; 
mu0=[0;0];

x0a=-0.1;
x0b=0.1;
wa=-0.01;
wb=0.01;
va=-0.01;
vb=0.01;

Sigmax0=eye(n)*1/3*(x0a^2+x0b^2+x0a*x0b);
Sigmaw=eye(n)*1/3*(wa^2+wb^2+wa*wb);
Sigmav=eye(p)*1/3*(va^2+vb^2+va*vb);



%%Stacked operator definition. "_b" stands for bold notation in the paper 
A_b=kron(eye(N+1),A);
B_b=[kron(eye(N),B);zeros(n,m*N)];
C_b=[kron(eye(N+1),C)];
Z=zeros(n*(N+1),n*(N+1));
for(i=1:N+1)
        for(j=1:N+1)
                if(i==j+1)
                        Z([(i-1)*n+1:i*n],[(j-1)*n+1:j*n])=eye(n);
                end
        end
end


P11=round(inv(eye(n*(N+1))-Z*A_b),3);
P12=round(inv(eye(n*(N+1))-Z*A_b)*Z*B_b,3);


Sigmaw_b=blkdiag(Sigmax0,kron(eye(N),Sigmaw));
Sigmav_b=kron(eye(N+1),Sigmav);
mu_w=[mu0;zeros(N*n,1)];



%%Cost function parameters
M=20*eye(p);
R=20*eye(m);
M_b=kron(eye(N+1),M);
R_b=kron(eye(N),R);



%INFO STRUCTURE

%struct = kron(tril(ones(N,(N+1))),[1 0 0]);
%struct=[0 0 0 0 0 0 0 0 0;1 0 0 0 0 0 0 0 0];
struct = [1 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0];


cardinality=sum(sum(struct));
K=sym('K',[m*(N),p*(N+1)]);
assume(K,'real');
K=K.*struct;


%%%stacking the non-zero decision variables into a single vector
vec_K=sym(zeros(cardinality,1));
count=0;
positions=zeros(cardinality,2);
for(i=1:m*N)
        for(j=1:p*(N+1))
                if(struct(i,j)==1)
                        count=count+1;
                        vec_K(count)=K(i,j);
                        positions(count,:)=[i,j];
                end
        end
end

OutputCost_convexprogram_outputfeedback;
OutputCost_create_cost;

ezsurf(cost) % to see that the cost is non-convex

grad = gradient(cost);

k=1;
gradient_dominance_poly = k*grad'*grad-cost+optimal_value; %% check that this is a SOS -> certified global gradient dominance with tau=1. I guess it is not well-behaved though... even if I reduce optimal_value it still finds a decomposition, even
[Q2,V2]=findsos(gradient_dominance_poly,'rational')