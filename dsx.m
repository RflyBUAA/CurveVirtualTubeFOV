function u=dsx(x,es)
A =-1/(32*es^5);
B =3/(16*es^5);
C =(5*(es^2 - 3))/(32*es^5);
D =-(5*(es - 1)*(es + 1))/(8*es^5);
E =-(15*(es - 1)^2*(es + 1)^2)/(32*es^5);
F =((es + 1)^3*(8*es^2 - 9*es + 3))/(16*es^5);
G =-((es - 1)^4*(5*es^2 + 4*es + 1))/(32*es^5);
 

if x<=1-es
    u = 1;
elseif 1-es<=x  && x<=1+es
    u = 6*A*x^5 + 5*B*x^4 + 4*C*x^3 + 3*D*x^2 + 2*E*x + F;
else
    u = 0;
end

end


 