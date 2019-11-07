%%rebuild K
Kiterate=zeros(m*N,p*N);
for(internal_i=1:cardinality)
        coord_x=positions(internal_i,1);
        coord_y=positions(internal_i,2);
        Kiterate(coord_x,coord_y)=parameters(internal_i);
end