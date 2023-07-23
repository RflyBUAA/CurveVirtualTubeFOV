function [vc,dvc,ddvc]  = CurveTubeVectorField(gamma,t,n,np,k,kr,Rgamma,rt,m,vmi,rs,ra,rd,alpha,lamba0,lambapi,dlamba0ds,dlambapids,dkds,dkrds,d2lamba0ds2,d2lambapids2,Pos,Vel,Acc)
global N
gamma_locate=zeros(1,N);
dti=zeros(1,N);
um1=zeros(2,N);
u2=zeros(2,N);
um3=zeros(2,N);
u=zeros(2,N);
umc=zeros(2,N);
vc=zeros(2,N);
cosmu=zeros(1,N);

eplision1 = 0.9;
eplision2 = 0.1;
k3 = 1.5;
k2 = 1.5;
%% Calculate Veloctiy Command
for i=1:N
    %% Calculate Nearest Point on the Generating Curve
    dis_gamma_min=10^8;
    for j=1:length(gamma(1,:))
        dis_gamma = norm(Pos(1:2,i)-gamma(:,j));
        if dis_gamma < dis_gamma_min
            gamma_locate(i) = j;
            dis_gamma_min = dis_gamma;
        end
    end
    %% Line Approaching Term
    um1(:,i) = vmi(i)*t(:,gamma_locate(i));
    %% Virtual Tube Keeping Term
    dti(i) = rt(gamma_locate(i))-norm(Pos(1:2,i)-m(:,gamma_locate(i)));
    um3(:,i) = -k3 * sigma2(dti(i),rs,ra) * (Pos(1:2,i)-m(:,gamma_locate(i)))/norm(Pos(1:2,i)-m(:,gamma_locate(i)));
    %% Robot Avoidance Term
    for j=1:N
        if i~=j
            pij = Pos(1:2,i) - Pos(1:2,j);
%             if norm(pij)<=rd && acos(dot(t(:,gamma_locate(i)),-pij)/(norm(t(:,gamma_locate(i)))*norm(-pij)))<alpha/2
            if norm(pij)<=rd && gamma_locate(i)<=gamma_locate(j)
                u2(:,i) = u2(:,i) + k2 * sigma2(norm(pij),2*rs,rs+ra)*(pij/norm(pij));
            end
        end
    end
    %% Calculation of u and umc
    u(:,i) = um1(:,i) + u2(:,i) + um3(:,i);
    umc(:,i) = sigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*u(:,i)/norm(u(:,i));
    %% Calculation of vc
    cosmu(i) = um1(:,i)'*umc(:,i)/(norm(um1(:,i))*norm(umc(:,i)));
    if dot(um1(:,i),umc(:,i))>=0
        vc(:,i) = (1-sigma2(cosmu(i),0,eplision2))*umc(:,i);
    else
        vc(:,i) = zeros(2,1);
    end
    if isnan(vc(:,i))
        vc(:,i)=zeros(2,1);
    end
end
%% Calculate First Derivative of Velocity Command with respect to Time
dum1=zeros(2,N);
du2=zeros(2,N);
dum3=zeros(2,N);
du=zeros(2,N);
dumc=zeros(2,N);
dvc=zeros(2,N);
ddti=zeros(1,N);
dcosmudt=zeros(1,N);
RgammaRtube=zeros(1,N);
Rtube=zeros(1,N);
for i=1:N
    %% Calculate Rgamma/Rtube
    Rtube(i) = Rgamma(gamma_locate(i)) - (Pos(1:2,i)-gamma(:,gamma_locate(i)))'*np(:,gamma_locate(i));
    RgammaRtube(i) = Rgamma(gamma_locate(i))/Rtube(i);
    if isnan(RgammaRtube(i))
        RgammaRtube(i) = 1;
    end
    %% Line Approaching Term
    dum1(:,i) = vmi(i)*kr(gamma_locate(i))*RgammaRtube(i)*t(:,gamma_locate(i))'*Vel(1:2,i)*n(:,gamma_locate(i));
    %% Virtual Tube Keeping Term
    drtdp = 1/2*(dlamba0ds(gamma_locate(i))+dlambapids(gamma_locate(i)))*RgammaRtube(i)*t(:,gamma_locate(i))';
    dmdp = ( t(:,gamma_locate(i)) + 1/2 * (dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i))) * n(:,gamma_locate(i))...
        -1/2 * (lamba0(gamma_locate(i))-lambapi(gamma_locate(i))) * kr(gamma_locate(i)) *t(:,gamma_locate(i)) )*RgammaRtube(i)*t(:,gamma_locate(i))';
    ddti(i) = (drtdp-(Pos(1:2,i)-m(:,gamma_locate(i)))'/norm(Pos(1:2,i)-m(:,gamma_locate(i)))*(eye(2)-dmdp))*Vel(1:2,i);

    dum3(:,i) = -k3*( dsigma2(dti(i),rs,ra) * ddti(i) * (Pos(1:2,i)-m(:,gamma_locate(i)))/norm(Pos(1:2,i)-m(:,gamma_locate(i)))...
        +sigma2(dti(i),rs,ra) * (Vel(1:2,i)-dmdp*Vel(1:2,i))/norm(Pos(1:2,i)-m(:,gamma_locate(i)))...
        -sigma2(dti(i),rs,ra) * (Pos(1:2,i)-m(:,gamma_locate(i)))*((Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i)))/((norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3)    );
    %% Robot Avoidance Term
    for j=1:N
        if i~=j
            pij = Pos(1:2,i) - Pos(1:2,j);
            vij = Vel(1:2,i) - Vel(1:2,j);
            if norm(pij)<=rd && acos(dot(t(:,gamma_locate(i)),-pij)/(norm(t(:,gamma_locate(i)))*norm(-pij)))<alpha/2
                du2(:,i) = du2(:,i) + k2*(dsigma2(norm(pij),2*rs,rs+ra)*pij'/norm(pij)*vij*pij/norm(pij)...
                    +sigma2(norm(pij),2*rs,rs+ra)*(vij*(norm(pij))^2-pij*pij'*vij)/((norm(pij))^3));
            end
        end
    end
    %% Calculation of du and dumc
    du(:,i) = dum1(:,i) + du2(:,i) + dum3(:,i);
    dumc(:,i) = dsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*u(:,i)'/norm(u(:,i))*du(:,i)*u(:,i)/norm(u(:,i))...
        +sigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*(du(:,i)*(norm(u(:,i)))^2-u(:,i)*u(:,i)'*du(:,i))/((norm(u(:,i)))^3);
    %% Calculation of dvc
    dcosmudt(i) = (dumc(:,i)'*um1(:,i)+umc(:,i)'*dum1(:,i))/(norm(umc(:,i))*norm(um1(:,i)))...
        -(umc(:,i)'*um1(:,i))*(umc(:,i)'*dumc(:,i)*(norm(um1(:,i)))^2+um1(:,i)'*dum1(:,i)*(norm(umc(:,i)))^2)/((norm(umc(:,i)))^3*(norm(um1(:,i)))^3);
    if dot(um1(:,i),umc(:,i))>=0
        dvc(:,i) = (1-sigma2(cosmu(i),0,eplision2))*dumc(:,i)-dsigma2(cosmu(i),0,eplision2)*dcosmudt(i)*umc(:,i);
    else
        dvc(:,i) = zeros(2,1);
    end
    if isnan(dvc(:,i))
        dvc(:,i)=zeros(2,1);
    end
end
%% Calculate Second Derivative of Velocity Command with respect to Time
ddum1=zeros(2,N);
ddu2=zeros(2,N);
ddum3=zeros(2,N);
ddu=zeros(2,N);
ddumc=zeros(2,N);
ddvc=zeros(2,N);
dRgammaRtube=zeros(1,N);
dddti=zeros(1,N);
d2cosmudt2=zeros(1,N);
for i=1:N
    %% Calculate Relative Variable
    dsi = RgammaRtube(i)*t(:,gamma_locate(i))'*Vel(1:2,i);
    dk = dkds(gamma_locate(i)) * dsi;
    dkr = dkrds(gamma_locate(i)) * dsi;
    dRgamma = -1/((k(gamma_locate(i)))^2) * dk;
    dRtube = dRgamma-(Vel(1:2,i)-t(:,gamma_locate(i))*dsi)'*np(:,gamma_locate(i))...
        +(Pos(1:2,i)-gamma(:,gamma_locate(i)))'*k(gamma_locate(i))*t(:,gamma_locate(i))*dsi;
    dRgammaRtube(i) = (Rtube(i)*dRgamma-dRtube*Rgamma(gamma_locate(i)))/((Rtube(i))^2);
    if isnan(dRgammaRtube(i))
        dRgammaRtube(i) = 0;
    end
    %% Line Approaching Term
    ddum1(:,i) = vmi(i)*dkr*RgammaRtube(i)*t(:,gamma_locate(i))'*Vel(1:2,i)*n(:,gamma_locate(i))...
        + vmi(i)*kr(gamma_locate(i))*dRgammaRtube(i)*t(:,gamma_locate(i))'*Vel(1:2,i)*n(:,gamma_locate(i))...
        + vmi(i)*(kr(gamma_locate(i)))^2*RgammaRtube(i)*dsi*n(:,gamma_locate(i))'*Vel(1:2,i)*n(:,gamma_locate(i))...
        + vmi(i)*kr(gamma_locate(i))*RgammaRtube(i)*t(:,gamma_locate(i))'*Acc(1:2,i)*n(:,gamma_locate(i))...
        - vmi(i)*(kr(gamma_locate(i)))^2*RgammaRtube(i)*dsi*t(:,gamma_locate(i))'*Vel(1:2,i)*t(:,gamma_locate(i));
    %% Robot Avoidance Term
    for j=1:N
        if i~=j
            pij = Pos(1:2,i) - Pos(1:2,j);
            vij = Vel(1:2,i) - Vel(1:2,j);
            aij = Acc(1:2,i) - Acc(1:2,j);
            if norm(pij)<=rd && acos(dot(t(:,gamma_locate(i)),-pij)/(norm(t(:,gamma_locate(i)))*norm(-pij)))<alpha/2
                ddu2(:,i) = ddu2(:,i)+k2*(ddsigma2(norm(pij),2*rs,rs+ra)*(pij'/norm(pij)*vij)^2*pij/norm(pij)...
                    + dsigma2(norm(pij),2*rs,rs+ra)*(vij'*(norm(pij))^2-pij'*(pij'*vij))/((norm(pij))^3)*vij*pij/norm(pij)...
                    + dsigma2(norm(pij),2*rs,rs+ra)*pij'/norm(pij)*aij*pij/norm(pij)...
                    + 2*dsigma2(norm(pij),2*rs,rs+ra)*pij'/norm(pij)*vij*(vij*(norm(pij))^2-pij*(pij'*vij))/((norm(pij))^3)...
                    + sigma2(norm(pij),2*rs,rs+ra)*(( (norm(pij))^2*(aij*(norm(pij))^2+vij*(pij'*vij)-pij*(vij'*vij)-pij*(pij'*aij))...
                    - 3*pij'*vij*(vij*(norm(pij))^2-pij*(pij'*vij)) )/((norm(pij))^5)) );
            end
        end
    end
    %% Virtual Tube Keeping Term
    drtdp = 1/2*(dlamba0ds(gamma_locate(i))+dlambapids(gamma_locate(i)))*RgammaRtube(i)*t(:,gamma_locate(i))';
    dmdp = ( t(:,gamma_locate(i)) + 1/2 * (dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i))) * n(:,gamma_locate(i))...
        -1/2 * (lamba0(gamma_locate(i))-lambapi(gamma_locate(i))) * kr(gamma_locate(i)) *t(:,gamma_locate(i)) )*RgammaRtube(i)*t(:,gamma_locate(i))';
    d2rtdpdt = 1/2*(d2lamba0ds2(gamma_locate(i))+d2lambapids2(gamma_locate(i)))*dsi*RgammaRtube(i)*t(:,gamma_locate(i))'...
        + 1/2*(dlamba0ds(gamma_locate(i))+dlambapids(gamma_locate(i)))*dRgammaRtube(i)*t(:,gamma_locate(i))'...
        + 1/2*(dlamba0ds(gamma_locate(i))+dlambapids(gamma_locate(i)))*RgammaRtube(i)*kr(gamma_locate(i))*n(:,gamma_locate(i))'*dsi;
    d2mdpdt = (kr(gamma_locate(i))*n(:,gamma_locate(i))*dsi ...
        + 1/2*(d2lamba0ds2(gamma_locate(i))-d2lambapids2(gamma_locate(i)))*dsi*n(:,gamma_locate(i))...
        - 1/2*(dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i)))*kr(gamma_locate(i))*t(:,gamma_locate(i))*dsi ...
        - 1/2*(dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i)))*dsi*kr(gamma_locate(i))*t(:,gamma_locate(i)) ...
        - 1/2*(lamba0(gamma_locate(i))-lambapi(gamma_locate(i)))*dkr*t(:,gamma_locate(i)) ...
        - 1/2*(lamba0(gamma_locate(i))-lambapi(gamma_locate(i)))*(kr(gamma_locate(i)))^2*n(:,gamma_locate(i))*dsi)*RgammaRtube(i)*t(:,gamma_locate(i))'...
        + (t(:,gamma_locate(i)) + 1/2 * (dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i))) * n(:,gamma_locate(i))...
        - 1/2 * (lamba0(gamma_locate(i))-lambapi(gamma_locate(i))) * kr(gamma_locate(i)) *t(:,gamma_locate(i)))*dRgammaRtube(i)*t(:,gamma_locate(i))'...
        + (t(:,gamma_locate(i)) + 1/2 * (dlamba0ds(gamma_locate(i))-dlambapids(gamma_locate(i))) * n(:,gamma_locate(i))...
        - 1/2 * (lamba0(gamma_locate(i))-lambapi(gamma_locate(i))) * kr(gamma_locate(i)) *t(:,gamma_locate(i)))*RgammaRtube(i)*kr(gamma_locate(i))*n(:,gamma_locate(i))'*dsi;
    dddti(i) = ( d2rtdpdt-(Vel(1:2,i)-dmdp*Vel(1:2,i))'/norm(Pos(1:2,i)-m(:,gamma_locate(i)))*(eye(2)-dmdp) ...
        +(Pos(1:2,i)-m(:,gamma_locate(i)))'*((Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i)))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3*(eye(2)-dmdp)...
        +(Pos(1:2,i)-m(:,gamma_locate(i)))'/norm(Pos(1:2,i)-m(:,gamma_locate(i)))*d2mdpdt )*Vel(1:2,i)...
        +(drtdp-(Pos(1:2,i)-m(:,gamma_locate(i)))'/norm(Pos(1:2,i)-m(:,gamma_locate(i)))*(eye(2)-dmdp))*Acc(1:2,i);
    ddum3(:,i) = -k3*(ddsigma2(dti(i),rs,ra)*(ddti(i))^2*(Pos(1:2,i)-m(:,gamma_locate(i)))/norm(Pos(1:2,i)-m(:,gamma_locate(i))) ...
        + dsigma2(dti(i),rs,ra)*dddti(i)*(Pos(1:2,i)-m(:,gamma_locate(i)))/norm(Pos(1:2,i)-m(:,gamma_locate(i)))...
        + dsigma2(dti(i),rs,ra)*ddti(i)*((Vel(1:2,i)-dmdp*Vel(1:2,i))/norm(Pos(1:2,i)-m(:,gamma_locate(i))) ...
        -(Pos(1:2,i)-m(:,gamma_locate(i)))*(Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3) ...
        + dsigma2(dti(i),rs,ra)*ddti(i)*(Vel(1:2,i)-dmdp*Vel(1:2,i))/norm(Pos(1:2,i)-m(:,gamma_locate(i))) ...
        + sigma2(dti(i),rs,ra)*(Acc(1:2,i)-d2mdpdt*Vel(1:2,i)-dmdp*Acc(1:2,i))/norm(Pos(1:2,i)-m(:,gamma_locate(i))) ...
        - sigma2(dti(i),rs,ra)*(Vel(1:2,i)-dmdp*Vel(1:2,i))*(Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3 ...
        - dsigma2(dti(i),rs,ra)*ddti(i)*(Pos(1:2,i)-m(:,gamma_locate(i)))*(Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3 ...
        - sigma2(dti(i),rs,ra)*(Vel(1:2,i)-dmdp*Vel(1:2,i))*(Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3 ...
        - sigma2(dti(i),rs,ra)*(Pos(1:2,i)-m(:,gamma_locate(i)))*(Vel(1:2,i)-dmdp*Vel(1:2,i))'*(Vel(1:2,i)-dmdp*Vel(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3 ...
        - sigma2(dti(i),rs,ra)*(Pos(1:2,i)-m(:,gamma_locate(i)))*(Pos(1:2,i)-m(:,gamma_locate(i)))'*(Acc(1:2,i)-d2mdpdt*Vel(1:2,i)-dmdp*Acc(1:2,i))/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^3 ...
        + sigma2(dti(i),rs,ra)*3*(Pos(1:2,i)-m(:,gamma_locate(i)))*((Pos(1:2,i)-m(:,gamma_locate(i)))'*(Vel(1:2,i)-dmdp*Vel(1:2,i)))^2/(norm(Pos(1:2,i)-m(:,gamma_locate(i))))^5 );
    %% Calculation of ddu and ddumc
    ddu(:,i) = ddum1(:,i) + ddu2(:,i) + ddum3(:,i);
    ddumc(:,i) = ddsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*(u(:,i)'/norm(u(:,i))*du(:,i))^2*u(:,i)/norm(u(:,i)) ...
        + dsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*(du(:,i)*(norm(u(:,i)))^2-u(:,i)*u(:,i)'*du(:,i))'/((norm(u(:,i)))^3)*du(:,i)*u(:,i)/norm(u(:,i)) ...
        + dsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*u(:,i)'/norm(u(:,i))*ddu(:,i)*u(:,i)/norm(u(:,i)) ...
        + dsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*u(:,i)'/norm(u(:,i))*du(:,i)*(du(:,i)*(norm(u(:,i)))^2-u(:,i)*u(:,i)'*du(:,i))/((norm(u(:,i)))^3) ...
        + dsigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*u(:,i)'/norm(u(:,i))*du(:,i)*(du(:,i)*(norm(u(:,i)))^2-u(:,i)*u(:,i)'*du(:,i))/((norm(u(:,i)))^3) ...
        + sigma1(norm(u(:,i)),eplision1*vmi(i),vmi(i))*((norm(u(:,i)))^2*(ddu(:,i)*(norm(u(:,i)))^2+du(:,i)*u(:,i)'*du(:,i)-u(:,i)*du(:,i)'*du(:,i)-u(:,i)*u(:,i)'*ddu(:,i)) ...
        - 3*u(:,i)'*du(:,i)*(du(:,i)*(norm(u(:,i)))^2-u(:,i)*u(:,i)'*du(:,i)))/((norm(u(:,i)))^5);
    %% Calculation of ddvc
    d2cosmudt2(i) = (ddumc(:,i)'*um1(:,i)+2*dumc(:,i)'*dum1(:,i)+umc(:,i)'*ddum1(:,i))/(norm(umc(:,i))*norm(um1(:,i)))...
        - 2*(dumc(:,i)'*um1(:,i)+umc(:,i)'*dum1(:,i))*(umc(:,i)'*dumc(:,i)*(norm(um1(:,i)))^2+um1(:,i)'*dum1(:,i)*(norm(umc(:,i)))^2)/((norm(umc(:,i)))^3*(norm(um1(:,i)))^3) ...
        - umc(:,i)'*um1(:,i)*(dumc(:,i)'*dumc(:,i)*(norm(um1(:,i)))^2+umc(:,i)'*ddumc(:,i)*(norm(um1(:,i)))^2+2*umc(:,i)'*dumc(:,i)*um1(:,i)'*dum1(:,i))/((norm(umc(:,i)))^3*(norm(um1(:,i)))^3) ...
        - umc(:,i)'*um1(:,i)*(dum1(:,i)'*dum1(:,i)*(norm(umc(:,i)))^2+um1(:,i)'*ddum1(:,i)*(norm(umc(:,i)))^2+2*um1(:,i)'*dum1(:,i)*umc(:,i)'*dumc(:,i))/((norm(umc(:,i)))^3*(norm(um1(:,i)))^3) ...
        + umc(:,i)'*um1(:,i)*(umc(:,i)'*dumc(:,i)*(norm(um1(:,i)))^2+um1(:,i)'*dum1(:,i)*(norm(umc(:,i)))^2)*(3*umc(:,i)'*dumc(:,i)*(norm(um1(:,i)))^2+3*um1(:,i)'*dum1(:,i)*(norm(umc(:,i)))^2)/((norm(umc(:,i)))^5*(norm(um1(:,i)))^5);
    if dot(um1(:,i),umc(:,i))>=0
        ddvc(:,i) = -dsigma2(cosmu(i),0,eplision2)*dcosmudt(i)*dumc(:,i) ...
            + (1-sigma2(cosmu(i),0,eplision2))*ddumc(:,i) ...
            - ddsigma2(cosmu(i),0,eplision2)*(dcosmudt(i))^2*umc(:,i) ...
            - dsigma2(cosmu(i),0,eplision2)*d2cosmudt2(i)*umc(:,i) ...
            - dsigma2(cosmu(i),0,eplision2)*dcosmudt(i)*dumc(:,i);
    else
        ddvc(:,i) = zeros(2,1);
    end
    if isnan(ddvc(:,i))
        ddvc(:,i)=zeros(2,1);
    end 
end

