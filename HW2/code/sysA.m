%% Build the system matrix
function [A,b] = sysA(rho, Lambda)
    beta = [0.060,0.364,0.349,0.628,0.179,0.070];
    lambda = [0.0129,0.0311,0.134,0.331,1.26,3.21];
   
    A = zeros(7,7);
    for i = 1:6
        A(i,i) = - lambda(i);
        A(i,7) = beta(i);
        A(7,i) = lambda(i)/Lambda(1);
    end
    A(7,7) = (rho(1)*sum(beta) - sum(beta) )/Lambda(1);
    
    b = [beta./lambda, 1];
end
