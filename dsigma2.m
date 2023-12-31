function u=dsigma2(x,d1,d2)
A =-20/(d1 - d2)^7; 
B =(70*(d1 + d2))/(d1 - d2)^7;
C =-(84*(d1^2 + 3*d1*d2 + d2^2))/(d1 - d2)^7;
D =(35*(d1^3 + 9*d1^2*d2 + 9*d1*d2^2 + d2^3))/(d1 - d2)^7;
E =-(140*d1*d2*(d1^2 + 3*d1*d2 + d2^2))/(d1 - d2)^7;
F =(210*d1^2*d2^2*(d1 + d2))/(d1 - d2)^7;
G =-(140*d1^3*d2^3)/(d1 - d2)^7;
H =(d2^4*(35*d1^3 - 21*d1^2*d2 + 7*d1*d2^2 - d2^3))/(d1 - d2)^7; 

if x<=d1
    u = 0;
elseif d1<=x  && x<=d2
    u = 7*A*x^6 + 6*B*x^5 + 5*C*x^4 + 4*D*x^3 + 3*E*x^2 + 2*F*x + G;
else
    u = 0;
end

end


 