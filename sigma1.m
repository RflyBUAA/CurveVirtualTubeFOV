function u=sigma1(x,d1,d2)
% A =3/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
% B =-(8*d1 + 7*d2)/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
% C =(2*(3*d1^2 + 10*d1*d2 + 2*d2^2))/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
% D =-(6*d1*(2*d2^2 + 3*d1*d2))/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
% E =(18*d1^2*d2^2 - 4*d1*d2^3 + d2^4)/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
% F =(d1^3*d2*(d1 - 4*d2))/((d1 - d2)*(d1^3 - 3*d1^2*d2 + 3*d1*d2^2 - d2^3));
A =3/(d1 - d2)^4;
B =-(8*d1 + 7*d2)/(d1 - d2)^4;
C =(6*d1^2 + 20*d1*d2 + 4*d2^2)/(d1 - d2)^4;
D =-(6*d1*d2*(3*d1 + 2*d2))/(d1 - d2)^4;
E =(d2^2*(18*d1^2 - 4*d1*d2 + d2^2))/(d1 - d2)^4;
F =(d1^3*d2*(d1 - 4*d2))/(d1 - d2)^4;

u=0;
if x<=d1
    u = x;
elseif x>d1 && x<d2
    u = A*x^5+B*x^4+C*x^3+D*x^2+E*x+F;
elseif x>=d2
    u = d2;
end
