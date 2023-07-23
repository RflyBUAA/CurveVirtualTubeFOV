function jc = VelocityController(vc,dvc,ddvc,Vel,Acc,Kv,Ka)
global N
jc=zeros(2,N);
for i=1:N
    jc(:,i)=Kv*(vc(:,i)-Vel(:,i))+Ka*(dvc(:,i)-Acc(:,i))+ddvc(:,i);
end
end