
%%Linear system definition
%A=[1.6 0 0 0 0;0.5 1.6 0 0 0;2.5 2.5 -1.4 0 0;-2 1 -2 0.1 0;0 2 0 -0.5 1.1];
A =  [1 0 -10;-1 1 0;0 0 1];


n=size(A,1);
B = [1;-1;0];

C = eye(n);

m=size(B,2);
p=size(C,1);

N=2; 
mu0=[0.1;-0.1;0.1];

x0a=-0.01;
x0b=0.01;
wa=-0.001;
wb=0.001;
va=-0.001;
vb=0.001;

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
M=1*eye(p);
R=1*eye(m);
M_b=kron(eye(N+1),M);
R_b=kron(eye(N),R);



%%Information structure definition
%S = [0 0 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 0;0 1 0 0 1];
%S=[1 0;1 1];
%S=1;
%struct_small=ones(N,N);
%struct_small=tril(struct_small);
%struct=kron(struct_small,S); %%%% this struct is \textbf{S} in the paper. QI can be verified by checking that the structure of struct*C_b*P12*struct is included in struct itself.

struct = kron(tril(ones(N,(N+1))),[1 0 0]);
%struct(end,:) = zeros(1,p*N); %because it does not make sense to have an input at the last time step, since it does not influence the last output.
%struct=zeros(m*N,p*N);
%struct(1,1)=1;

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

