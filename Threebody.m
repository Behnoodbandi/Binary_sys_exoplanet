function [X,V,Kinetic,Potential,Total_energy]=Threebody(W,M,mp,G,X0,V_p,dt,R,Max_step,dx)
%Test Inputs W=G=R=1 M=1000 mp=0.001  
Output=zeros(Max_step,9);
figure;
hold;
X00=X0+rand(1,3)*dx;
delta0=sqrt((X0-X00)*((X0-X00)'));
plot3(R*cos(0:0.01:2*pi),R*sin(0:0.01:2*pi),zeros(size(0:0.01:2*pi,2)));
X=zeros(Max_step,3);
X1=zeros(Max_step,3);
V=zeros(Max_step,3);
Kinetic=zeros(Max_step-1,1);
Potential=zeros(Max_step-1,1);
Total_energy=zeros(Max_step-1,1);
X(1,:)=X0;
X1(1,:)=X00;
ax =gca;
title(['X0=' , '[' num2str(X(1,1)) ',' num2str(X(1,2)) ','  num2str(X(1,3)) '] V0=[' num2str(V_p(1)) ',' num2str(V_p(2)) ',' num2str(V_p(3)) ']' ',' num2str(dt)])
Name=['X0=' , '[' num2str(X(1,1)) ',' num2str(X(1,2)) ','  num2str(X(1,3)) '] V0=[' num2str(V_p(1)) ',' num2str(V_p(2)) ',' num2str(V_p(3)) ']'];
for i=1:Max_step-1
    r1=R*[cos(W*(i-1)*dt),sin(W*(i-1)*dt),0];
    r2=-r1;
    N1=sqrt((X(i,:)-r1)*(X(i,:)-r1)');
    N2=sqrt((X(i,:)-r2)*(X(i,:)-r2)');
    Kinetic(i)=(mp*(V_p*V_p'))/2;
    Potential(i)=-(G*M*mp)*((1/N1)+(1/N2));
    Total_energy(i)=Kinetic(i)+Potential(i);
    a_p=-(G*M)*((X(i,:)-r1)/(N1^3)+(X(i,:)-r2)/(N2^3));
    V_p=V_p+dt*a_p;
    V(i,:)=V_p;
    X(i+1,:)=X(i,:)+dt*V_p;
    X1(i+1,:)=X1(i,:)+dt*V_p;
    if rem(i,100000)==0
        i
    end
end
plot3(X(:,1),X(:,2),X(:,3));
%plot3(X1(:,1),X1(:,2),X1(:,3));
R=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
grid on;
Lim=[-20,20];
ax.YLim=Lim;
ax.XLim=Lim;
ax.ZLim=Lim;
figure;
hold;
title(['Energy Graph , X0=' , '[' num2str(X(1,1)) ',' num2str(X(1,2)) ',' num2str(X(1,3)) ']' ])
grid on;
plot(1:Max_step-1,Total_energy)
plot(1:Max_step-1,Kinetic)
plot(1:Max_step-1,Potential)
%lyapanov 
figure
hold;
title(['Lyapabov Exp. , X0=' , '[' num2str(X(1,1)) ',' num2str(X(1,2)) ',' ',' num2str(X(1,3)) ']' ])
XX(:)=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
XX1(:)=sqrt(X1(:,1).^2+X1(:,2).^2+X1(:,3).^2);
plot(1:Max_step,log((XX-XX1)/delta0));
Output(:,1:3)=X;
Output(:,4:6)=V;
Output(1:end-1,7)=Kinetic;
Output(1:end-1,8)=Potential;
Output(1:end-1,9)=Total_energy;
csvwrite([Name '.csv'],Output)
end