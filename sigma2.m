function u=sigma2(x,d1,d2)
A =-20/(d1 - d2)^7; 
B =(70*(d1 + d2))/(d1 - d2)^7;
C =-(84*(d1^2 + 3*d1*d2 + d2^2))/(d1 - d2)^7;
D =(35*(d1^3 + 9*d1^2*d2 + 9*d1*d2^2 + d2^3))/(d1 - d2)^7;
E =-(140*d1*d2*(d1^2 + 3*d1*d2 + d2^2))/(d1 - d2)^7;
F =(210*d1^2*d2^2*(d1 + d2))/(d1 - d2)^7;
G =-(140*d1^3*d2^3)/(d1 - d2)^7;
H =(d2^4*(35*d1^3 - 21*d1^2*d2 + 7*d1*d2^2 - d2^3))/(d1 - d2)^7; 

if x<=d1
    u = 1;
elseif d1<=x  && x<=d2
    u = A*x^7+B*x^6+C*x^5+D*x^4+E*x^3+F*x^2+G*x+H ;
else
    u = 0;
end

end


 