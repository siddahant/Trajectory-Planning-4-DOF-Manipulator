clear; clc;
syms a_1 a_2 d3 t1 t2 t4 d_4 
a1 = 0.5;
a2 = 0.3;
d4 = 0.2;
% scara robot DH Table theta d a alpha
dh_scara_type=[t1, 0, a_1, 0;...
              t2, 0, a_2, 0;...
              0, d_4+d3, 0, 0;...
              t4, 0, 0, 0];

% End effactor matrix
T=simplify(DH(dh_scara_type));
xe=[T([1:3],4);[t1+t2+t4]]; % [PX PY PZ PSI=[t1+t2+t4]]
analytica_jacobian=[diff(xe,t1) diff(xe,t2) diff(xe,d3) diff(xe,t4)];

syms t
% Desired end-effector pose trajectory
p_d = [0.4*cos(t*pi/10); 0.4*sin(t*pi/10); 0.1*(1+sin(t))]; % positon
phi_d = t*pi/10 + 5*pi/12; % orientation

% Derivatives of the desired end-effector pose trajectory
p_ddot = diff(p_d);
phi_ddot = diff(phi_d);

delta_t = 0.01;
te = 5;
time=0:delta_t:te;

p_d=double(subs(p_d,t,time));
phi_d=double(subs(phi_d,t,time));
p_ddot = double(subs(p_ddot,t,time));
phi_ddot = double(subs(phi_ddot,t,time));

% Initial condition
q(:,1) = [pi/4, pi/4, 0.3, 0]; % t1 t2 d3 psi=[t1+t2+t3]
%% Algorithm
K=15*eye(4); %  positive definite (usually diagonal)  Gain matrix
for i = 1:length(time)    
 x_e(:,i)= double(subs(xe,{t1 t2 d3 t4 a_1 a_2 d_4 },{q(1,i) q(2,i) q(3,i) q(4,i) a1 a2 d4}));
 J_A = double(subs(analytica_jacobian,{t1 t2 a_1 a_2 },{q(1,i) q(2,i) a1 a2})); % Jacobian
 desired_position = [p_d(:,i); phi_d(i)];
 e(:,i) = desired_position - x_e(:,i);  % operational space error
 x_ddot=[p_ddot(:,i); phi_ddot(i)];
 q_dot = inv(J_A)*(x_ddot+K*e(:,i)); % Jacobian (Pseudo-)inverse
 q(:,i+1)=q(:,i)+delta_t*q_dot; % add new coloum q 
end

%% Plots comparing the desired end-effector pose with the actual
figure()
plot(time,p_d(1,:),'--');
hold on;
plot(time,x_e(1,:));
legend('Desired','Actual');
xlabel('Time in sec.');
ylabel('Px');
title('Position of x');
hold off
figure()
plot(time,p_d(2,:),'--');
hold on;
plot(time,x_e(2,:));
ylabel('Py');
xlabel('Time in sec.'); 
title('Position y');
hold off
figure()
plot(time,p_d(3,:),'--');
hold on;
plot(time,x_e(3,:));
title('Position z:');
xlabel('Time in sec'); 
ylabel('Pz');
hold off
figure();

plot(time,phi_d(1,:),'--');
hold on;
plot(time,x_e(4,:));
xlabel('Time in sec.');
ylabel('Phi in radian');
title('Phi')
hold off

function dh = DH(DH_Parameter_Matrix)
syms  theta d a alpha1
A = [cos(theta),-cos(alpha1)*sin(theta),sin(theta)*sin(alpha1),a*cos(theta);
     sin(theta),cos(theta)*cos(alpha1),-cos(theta)*sin(alpha1),a*sin(theta);
     0         , sin(alpha1)          , cos(alpha1)           ,d           ;
     0         , 0                   , 0                    ,1           ];

%%
n=size(DH_Parameter_Matrix,1);
d0=DH_Parameter_Matrix; 
    
dh=eye(4);

for i=1:n
  dh_next= subs(A,{theta d a alpha1},{d0(i,1),d0(i,2),d0(i,3),d0(i,4)});
  dh=dh*dh_next;% some equation
end
end