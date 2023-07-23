clear;
clc;
%% Robot Init
global N gamma t n np k kr Rgamma rt m vmi rs ra rd alpha lamba0 lambapi dlamba0ds dlambapids dkds dkrds d2lamba0ds2 d2lambapids2 Kv Ka

N=20; % Number
vmi=2.1*ones(N,1); % Maximum Speed
rs=0.2; % Safety Radius
ra=0.7; % Avoidance Radius
rd=2;   % Detection Radius
alpha=1.2*pi; % FOV

UAV_PInit=zeros(N,2); % Robot Initial Position
UAV_PInit(1,:)=[2,0];UAV_PInit(2,:)=[2,1];UAV_PInit(3,:)=[2,2];UAV_PInit(4,:)=[2,3];UAV_PInit(5,:)=[2,4];
UAV_PInit(6,:)=[3,0];UAV_PInit(7,:)=[3,1];UAV_PInit(8,:)=[3,2];UAV_PInit(9,:)=[3,3];UAV_PInit(10,:)=[3,4];
UAV_PInit(11,:)=[4,0];UAV_PInit(12,:)=[4,1];UAV_PInit(13,:)=[4,2];UAV_PInit(14,:)=[4,3];UAV_PInit(15,:)=[4,4];
UAV_PInit(16,:)=[5,0];UAV_PInit(17,:)=[5,1];UAV_PInit(18,:)=[5,2];UAV_PInit(19,:)=[5,3];UAV_PInit(20,:)=[5,4];
UAV_PInit=UAV_PInit+[0.4,0.2];
UAV_VInit=zeros(N,2); % Robot Initial Velocity
UAV_AInit=zeros(N,2); % Robot Initial Acceleration

UAV_Initial=[UAV_PInit(:,1);UAV_PInit(:,2);UAV_VInit(:,1);UAV_VInit(:,2);UAV_AInit(:,1);UAV_AInit(:,2)];
%% Velocity  Controller Parameter
Kv = 5;
Ka = 5;

%% Generating Curve
curveA=2;
curveB=0.1*pi;
curve_x=0:0.04:30;
curve_y=curveA*sin(curveB*curve_x);
gamma=[curve_x;curve_y]; % Generating Curve

t=zeros(2,length(gamma(1,:)));
n=zeros(2,length(gamma(1,:)));
np=zeros(2,length(gamma(1,:)));
kr=zeros(1,length(gamma(1,:)));
k=zeros(1,length(gamma(1,:)));
Rgamma=zeros(1,length(gamma(1,:)));
dcurve_x=zeros(1,length(gamma(1,:)));
dcurve_y=zeros(1,length(gamma(1,:)));
ddcurve_x=zeros(1,length(gamma(1,:)));
ddcurve_y=zeros(1,length(gamma(1,:)));
for i=1:length(gamma(1,:))
    dcurve_x(i)=1;
    dcurve_y(i)=curveA*curveB*cos(curveB*(curve_x(i)));
    t(:,i)=[dcurve_x(i);dcurve_y(i)]/norm([dcurve_x(i);dcurve_y(i)]); % Unit Tangent Vector
    n(:,i)=[-t(2,i);t(1,i)]; % Unit Normal Vector
    ddcurve_x(i)=0;
    ddcurve_y(i)=-curveA*curveB^2*sin(curveB*(curve_x(i)));
    kr(i)=(dcurve_x(i)*ddcurve_y(i)-ddcurve_x(i)*dcurve_y(i))/(sqrt( (dcurve_x(i)^2+dcurve_y(i)^2)^3 ) ); % Relative Curvature
    k(i)=abs(kr(i)); % Curvature
    Rgamma(i)=1/k(i); % Radius of Curvature
    np(:,i)=(kr(i)/k(i))*n(:,i); % Principle Unit Normal Vector
end

%% Curve Virtual Tube
lamba0A=0.8;
lamba0B=0.15*pi;
lamba0C=2.5;
lamba0=lamba0A*sin(lamba0B*curve_x)+lamba0C; % Left Radius
lambapiA=0.8;
lambapiB=0.2*pi;
lambapiC=3;
lambapi=lambapiA*sin(lambapiB*curve_x)+lambapiC; % Right Radius

boundary_left=zeros(2,length(gamma(1,:)));
boundary_right=zeros(2,length(gamma(1,:)));
for i=1:length(gamma(1,:))
    boundary_left(:,i)=gamma(:,i)+lamba0(i)*n(:,i); % Left Boundary
    boundary_right(:,i)=gamma(:,i)-lambapi(i)*n(:,i); % Right Boundary
end

rt=0.5*(lamba0+lambapi); % Width
m=zeros(2,length(gamma(1,:))); % Middle Point
for i=1:length(gamma(1,:))
    m(:,i)=gamma(:,i)+0.5*(lamba0(i)-lambapi(i))*n(:,i);
end
%% High Order Calculation
dlamba0ds=zeros(1,length(gamma(1,:))); % \frac{\partial \lambda\left(s,0\right)}{\partial s}
dlambapids=zeros(1,length(gamma(1,:))); % \frac{\partial \lambda\left(s,\pi\right)}{\partial s}
dkds=zeros(1,length(gamma(1,:))); % \frac{\partial k}{\partial s}
dkrds=zeros(1,length(gamma(1,:))); % \frac{\partial k_r}{\partial s}
for i=1:length(gamma(1,:))-1
    ds=norm(gamma(:,i+1)-gamma(:,i));
    dlamba0ds(i)=(lamba0(i+1)-lamba0(i))/ds;
    dlambapids(i)=(lambapi(i+1)-lambapi(i))/ds;
    dkds(i)=(k(i+1)-k(i))/ds;
    dkrds(i)=(kr(i+1)-kr(i))/ds;
end
dlamba0ds(end)=dlamba0ds(end-1);
dlambapids(end)=dlambapids(end-1);
dkds(end)=dkds(end-1);
dkrds(end)=dkrds(end-1);

d2lamba0ds2=zeros(1,length(gamma(1,:))); % \frac{\partial^2 \lambda\left(s,0\right)}{\partial s^2}
d2lambapids2=zeros(1,length(gamma(1,:))); % \frac{\partial^2 \lambda\left(s,\pi\right)}{\partial s^2}
for i=1:length(gamma(1,:))-1
    ds=norm(gamma(:,i+1)-gamma(:,i));
    d2lamba0ds2(i)=(dlamba0ds(i+1)-dlamba0ds(i))/ds;
    d2lambapids2(i)=(dlambapids(i+1)-dlambapids(i))/ds;
end
d2lamba0ds2(end)=d2lamba0ds2(end-1);
d2lambapids2(end)=d2lambapids2(end-1);

%% Simulation
T=[0 12]; % Simulation Time
[tt,data] = ode45(@(tt,y) OdeThreeIntegrator(tt,y), T, UAV_Initial);

%% Result Plot
for i=1:20:length(data)
    hold on;
    axis([-5,45,-10,10])
    axis equal;
    plot(gamma(1,:),gamma(2,:),'k--','Linewidth',0.5);
    plot(boundary_left(1,:),boundary_left(2,:),'k','Linewidth',1);
    plot(boundary_right(1,:),boundary_right(2,:),'k','Linewidth',1);
    for j=1:N
        pos=[data(i,j);data(i,N+j)];
        mycircle_r(pos,rs);
    end
    pause(0.005);
    hold off
    plot(0,0);
end