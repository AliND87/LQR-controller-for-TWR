%Optimal Control Project
%Ali Nosouhi Dehnavi
%matricola:1950716
%Two wheeled robot LQR problem

clear all
close all
t=0:0.05:10;
%initial state
xint=[20; % initial pitch angle , 20 degree 
      0]; %  initial picth rate  , zero
  m=0.1; %kg Mass of the robot body
  l=20;% cm  Length of the robot body
  I=13.408;%kg.cm2 Inertia of the robot body
  g=9.806% m/s^2
  
%system dynamics
A=[0     1;
   (m*g*l)/(2*I)  0]; %0.73
B=[  0  ; 
   1/I];%0.074
C=[1 0; %to obtain second output change to [0 1]
  0 1];
D=[0;0];
sys1=ss(A,B,C,D);
Eig1=eig(A); %[0.85;-0.85] unstable
Co = ctrb(sys1)
length(A) - rank(Co) %== 0 system is controllable
Ob = obsv(sys1)
length(A) - rank(Ob)  % == 0 system is observable

%% To use in simulink
Q1=[1  0;  
   0  1];
R1=1;
[K1,S1,CLP1] = lqr(A,B,Q1,R1);
sys1_lqr=ss((A-B*K1),B,C,D);  %%%B  [0;0]
[y1,t,x1]=initial(sys1_lqr,xint,t);

%observer design
C1=[0 1];  % x2=angular velocity is the only available state %[0 0;0 1];
L=place(A',C1',[-8 -9])'; %Luenberger gain matrix , arbitrary values -8,-9
eig(A-L*C1)

Aobs=A-L*C1;
Bobs=[B L];
Cobs=C;% we want all estimated states
Dobs=zeros(2);
%% 



sw =input('case1(FULL STATE) enter 1 , case2(ESTIMATED STATE) enter 2 : ');
switch sw
    
 %Case1:LQR optimization for full state
 case 1

n=input('enter the number of optimization case = ');
Q=zeros(2,2,n);
R=zeros(1,n);
K=zeros(1,2,n);
S=zeros(2,2,n);
CLP=zeros(2,1,n);
% y=zeros(length(t),n);
% x=zeros(length(t),2,n);
% sys2=zeros(n);
for i=1:n

W=input(['Optimization parameters   '  num2str(i) ':',...
    '\nif you want to enter value of Q manually enter 1 else 2 = ']);
if W==1
    Q(1,1,i)=input('Please enter the weight value for angular(Q11) = ');
    Q(2,2,i)=input('Please enter the weight value for angular velocity(Q22) = ');
    Q(1,2,i)=0;
    Q(2,1,i)=0;
else
    Q(:,:,i)=transpose(C)*C;
end

R(i)=input('Please enter the weight value for control-torque R= ');
[K(:,:,i),S(:,:,i),CLP(:,:,i)] = lqr(A,B,Q(:,:,i),R(i));
 sys2(:,:,i)=ss((A-B*K(:,:,i)),[0;0],C,D);  %ss((A-B*K(:,:,i)),B,C,D)
  [y(:,:,i),t,x(:,:,i)]=initial(sys2(:,:,i),xint,t);
   u(:,i)=-(1/100)*K(:,:,i)*x(:,:,i)';
end

figure(1)
hold on
 
for r = 1 : n
    
    plot(t,y(:,1,r), '--', 'MarkerSize', 10, 'LineWidth', 2);
    legends{r} = sprintf(['Contr#%d Q=' num2str(Q(1,1,r)) ',' num2str(Q(2,2,r)) ' R=' num2str(R(r))], r);
    
end

legend(legends) % Display all the legend texts.
xlabel('$\textit{time [s]}$','Interpreter','latex','FontSize',14)
title('pitching angle [degree]','FontName','courier',...
      'FontSize',14)
grid on;

% 
figure(2)
hold on
 
for p = 1 : n
    
    plot(t,y(:,2,p), '--', 'MarkerSize', 10, 'LineWidth', 2);
    legends{p} = sprintf(['Contr#%d Q=' num2str(Q(1,1,p)) ',' num2str(Q(2,2,p)) ' R=' num2str(R(p))], p);
end
legend(legends) % Display all the legend texts.
xlabel('$\textit{time [s]}$','Interpreter','latex','FontSize',14)
title('pitching rate [deg/s]','FontName','courier',...
     'FontSize',14)
grid on;
% % 
figure(3)
hold on
 
for b = 1 : n
    
    plot(t,u(:,b), '--', 'MarkerSize', 10, 'LineWidth', 2);
    legends{b} = sprintf(['Contr#%d Q=' num2str(Q(1,1,b)) ',' num2str(Q(2,2,b)) ' R=' num2str(R(b))], b);
end
legend(legends) % Display all the legend texts.
xlabel('$\textit{time [s]}$','Interpreter','latex','FontSize',14)
title('torque [Nm]','FontName','courier',...
     'FontSize',14)
grid on;
% 






%% LQR optimization for stimated states

    case 2 
%% Luenberger Observer

%simple LQR wieght matrices for evaluate the effectiveness of 
%designed observer 
% case 2




 %  constructing new system composed of primary system and observer
% % % [  x1   ] 
% % % |  x2   |
% % % |x1(hat)|  = A_new . X + B_new .U
% % % [x2(hat)]

  A_new=[A -B*K1;
         L*C1 A-L*C1-B*K1];

  B_new=[0;0;0;0];  %%%[B;B];

     
  C_new=[C zeros(2);
         zeros(2) C ];

   D_new = [0;0;0;0]; 



 % closed loop system
  sys2=ss(A_new,B_new,C_new,D_new)
  eig(A_new)
  [y2,t,x2]=initial(sys2,[20;0;0;0],t);


  W1=input('if you want to enter value of Q manually enter 1 else 2 = ');
if W1==1
    Q3(1,1)=input('enter value of q11 = ');
    Q3(2,2)=input('enter value of q22 = ');
    Q3(1,2)=0;
    Q3(2,1)=0;
else
    Q3=transpose(C)*C;
end
% 
R3=input('enter the matrix R= ');
[K3,S3,CLP3] = lqr(A,B,Q3,R3);
 sys3=ss((A-B*K3),[0;0],C,D);   %ss((A-B*K3),B,C,D
 [y3,t,x3]=initial(sys3,xint,t);
  u3=(-K3*x3')/100;
  %using extended dynamic with observer
 A_ext=[A -B*K3;
         L*C1 A-L*C1-B*K3];  %A-L*C
 sys4=ss(A_ext,B_new,C_new,D_new)
 [y4,t,x4]=initial(sys4,[xint;0;0],t);
  u4=(-K3*x4(:,3:4)')/100;
    
  
figure(5)
Plot10=plot(t,y3(:,1),...
    t,y4(:,1)); grid
set(Plot10(1:2),'LineWidth',3);
yline(0,'-k');
xlabel('t(sec)','Interpreter','latex','FontSize',14)
title('full state vs observer pitch angle','FontName','courier','FontSize',14)
legend(Plot10,'x1-full','x1-observer','FontName','courier','FontSize',14)

figure(6)
Plot6=plot(t,y3(:,2),...
    t,y4(:,2)); grid
set(Plot6(1:2),'LineWidth',3);
yline(0,'-k');
xlabel('t(sec)','Interpreter','latex','FontSize',14)
title('full state vs observer pitch angle rate','FontName','courier','FontSize',14)
legend(Plot6,'x2-full','x2-observer','FontName','courier','FontSize',14)

figure(7)
Plot7=plot(t,u3,...
    t,u4); grid
set(Plot7(1:2),'LineWidth',3);
yline(0,'-k');
xlabel('t(sec)','Interpreter','latex','FontSize',14)
title('full state vs observer control effort','FontName','courier','FontSize',14)
legend(Plot7,'u-full','u-observer','FontName','courier','FontSize',14)


end

 
%  %% zero steady state error
%  dcgain(tf(sys1_lqr))
%  V=-inv((C*inv(A-B*K1))*B);% Gain matrix to assure y/r=1()when s-->0)
% 


