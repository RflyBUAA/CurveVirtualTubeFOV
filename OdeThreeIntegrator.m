function dydt = OdeThreeIntegrator(tt,y)
    global N gamma t n np k kr Rgamma rt m vmi rs ra rd alpha lamba0 lambapi dlamba0ds dlambapids dkds dkrds d2lamba0ds2 d2lambapids2 Kv Ka

    Pos =[y(1:N)';y(N+1:2*N)'];
    Vel = [y(2*N+1:3*N)';y(3*N+1:4*N)'];
    Acc = [y(4*N+1:5*N)';y(5*N+1:6*N)'];
    [vc,dvc,ddvc]  = CurveTubeVectorField(gamma,t,n,np,k,kr,Rgamma,rt,m,vmi,rs,ra,rd,alpha,lamba0,lambapi,dlamba0ds,dlambapids,dkds,dkrds,d2lamba0ds2,d2lambapids2,Pos,Vel,Acc);
    jc = VelocityController(vc,dvc,ddvc,Vel,Acc,Kv,Ka);
    u=[jc(1,:)';jc(2,:)'];
    dydt = zeros(6*N,1);
    dydt(1:4*N) = y(2*N+1:6*N);
    dydt(4*N+1:6*N) = u;
end