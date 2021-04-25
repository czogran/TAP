

M=[ 0 0;
    2 0;
    3 2;
    3.5 3;
    4 3.5];

psi= eye(5);

lambda = [ 0.25 0;
          0 0.25];
      
K=(M' * psi * M + lambda)^(-1)*M'