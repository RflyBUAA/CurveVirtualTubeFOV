function u=mysat(x,a)
% ensure the direction to be the same
    if norm(x)>a
       u =a*x/norm(x);
    else
       u =x ;
    end
    