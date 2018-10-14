clear
clc

global m1 m2 I1 I2 l1 l2 lg1 lg2 g Omega V

% �U�q1�͓V��ɒ����Ă���U�q(��͌��ʂŐ�)�A�U�q2�͐U�q1�̉��ɂ��Ă���U�q(��͌��ʂŐ�)

m1=385e-3;% �U�q1�̎���
m2=300e-3;% �U�q2�̎���
I1=5.54e-4;% �U�q1�̊������[�����g
I2=4.00e-4;% �U�q2�̊������[�����g
lg1=46e-3;% �U�q1�̉�]�΋􂩂�d�S�܂ł̋���
lg2=44e-3;% �U�q2�̉�]�΋􂩂�d�S�܂ł̋���
l1=83.4e-3;% �U�q1�̑S��
l2=100e-3;% �U�q2�̑S��

dt=1e-3; % �T���v�����O����
tmax=10;
kmax=tmax/dt;% �T���v���_��


% �����p�ψ�(deg)�C�������C�f�[�^�o�͂�rad
theta0_1=0;
theta0_2=0;

% �����p���x(rad/s)
dtheta0_1=0;
dtheta0_2=0;

Omega=2*2*pi; % ���U���g��
V=0.002;      % ��d��

kakudo=-pi/2; % ��Ίp�̉^���������ł��邽��(������0rad)
theta0_1=theta0_1/(180/pi)+kakudo;
theta0_2=theta0_2/(180/pi);
theta0_2=theta0_1+theta0_2;

g=9.8; % �d�͉����x

% ��͂��߂ɕK�v�Ȃ���
x=zeros(4,1);
x=[dtheta0_1;dtheta0_2;theta0_1;theta0_2];
logx=zeros(kmax,4);
logxs=logx;

% �����Q�N�b�^�@�ŋߎ�����
for k=1:kmax
    logx(k,:)=x';
    t=(k-1)*dt;
    k1=EoMInvPendForced2dof(x,t);
    k2=EoMInvPendForced2dof(x+k1*dt/2,t+dt/2);
    k3=EoMInvPendForced2dof(x+k2*dt/2,t+dt/2);
    k4=EoMInvPendForced2dof(x+k3*dt,t+dt);
    x=x+(k1+2*k2+2*k3+k4)/6*dt;    
end

t=(0:kmax-1)'*dt;

for k=1:kmax
    if mod(k,20)==1
        x=logx(k,:)';
        xs=logxs(k,:)';
        figure(2)
        plot([0 l1*cos(x(3))],[0 l1*sin(x(3))],'bo-',...
            [l1*cos(x(3)) l1*cos(x(3))+l2*cos(x(4))],[l1*sin(x(3)) l1*sin(x(3))+l2*sin(x(4))],'ro-');
        axis equal;
        xlim([-0.2 0.2]);
        ylim([-0.2 0.2]);
       drawnow; 
    end
end

index=(kmax-2/dt+1):kmax;
theta1=logx(index,3)-kakudo;
theta2=(logx(index,4)-logx(index,3));
dlmwrite('calc_data.csv',[t(index),theta1,theta2],'delimiter','\t','precision','%17.15e');
figure(3)
plot(t(index),[(theta1*180)./pi (theta2*180)./pi]);
xlabel('���� [s]');
ylabel('�p�ψ� [deg]');
legend('\theta_1','\theta_2');