clc,clear,close all

%% Cerinta 1

theta1_init=0.1;
theta2_init=0.2;

%% Cerinta 2

% modelul simulink

%% Cerinta 3
% in raport

%% Cerinta 4

% Parametrii sistemului
M=1; L=1; g=9.81; k=100;
F=0;

% Funcția f(t, x)
f = @(t, x) [
    x(3); % d_x1 = d_theta1 = x3
    x(4); % d_x2 = d_theta2 = x4
    (1/(M*L^2))*(F-M*g*L*sin(x(1))-(1/4)*k*L^2*(sin(x(1))-sin(x(2)))*cos(x(1))); % d_x3 = dd_theta1
    (1/(M*L^2))*(-M*g*L*sin(x(2))+(1/4)*k*L^2*(sin(x(1))-sin(x(2)))*cos(x(2))) % d_x4 = dd_theta2
];

h = 0.1; 
t = 0:h:60; 
x0 = [theta1_init; theta2_init; 0; 0]; % Condiții inițiale [theta1,  theta2, theta1_dot, theta2_dot]

x = zeros(4, length(t));
x(:, 1) = x0;

% Metoda Runge-Kutta
for i = 1:(length(t)-1)
    k1 = f(t(i), x(:, i));
    k2 = f(t(i) + h/2, x(:, i) + (h/2) * k1);
    k3 = f(t(i) + h/2, x(:, i) + (h/2) * k2);
    k4 = f(t(i) + h, x(:, i) + h * k3);
    x(:, i+1) = x(:, i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

figure;
plot(t, x(1, :), 'b', t, x(2, :), 'r');
xlabel('Timp (s)');
ylabel('Theta (rad)');
legend('\theta_1', '\theta_2');
title('Evoluția unghiurilor \theta_1 și \theta_2');

%% Cerinta 5 

figure
set(gcf,'Position',[100, 100, 1200, 300])
subplot(1,3,1)
plot(x(1,:),x(3,:))
title("dtheta1 in functie de theta1")
xlabel("theta1 (rad)")
ylabel("dtheta1 (rad/s)")

subplot(1,3,2)
plot(x(2,:),x(4,:))
title("dtheta2 in functie de theta2")
xlabel("theta2 (rad)")
ylabel("dtheta2 (rad/s)")

subplot(1,3,3)
plot(x(1,:),x(2,:))
title("theta2 in functie de theta1")
xlabel("theta1 (rad)")
ylabel("theta2 (rad)")

%% Cerinta 6 

y_RK=x(1:2,:);

Tmax=60;
u=zeros(size(t));
F=timeseries(u,t);

load_system("model_simulink")
set_param("model_simulink","StopTime",num2str(Tmax))
out=sim("model_simulink");
y_Slx=[out.theta1.Data'; out.theta2.Data'];

t_Slx=out.tout;
y_RK_interp = interp1(t, y_RK', t_Slx);

y_RK_interp=y_RK_interp';

diferenta=(y_RK_interp - y_Slx)';
eroare = norm(y_RK_interp - y_Slx); 
disp(['Eroarea de integrare (norma 2): ', num2str(eroare)]);

figure;
plot(t_Slx, diferenta(:, 1), 'b', t_Slx, diferenta(:, 2), 'r');
xlabel('Timp (s)');
ylabel('Diferență');
legend('Diferență theta_1', 'Diferență theta_2');
title('Eroarea de integrare între Simulink și Runge-Kutta');
grid on;

%% Cerinta 7 

i=1;
for k=1:0.1:5
    u=k.*double(t>=0);
    F=timeseries(u,t);
    out=sim("model_simulink");
    ystar(i,1)=out.theta1.Data(end);
    ystar(i,2)=out.theta2.Data(end);
    Fstar(i)=k;
    i=i+1;
end
plot(Fstar,ystar(:,1),'x',Fstar,ystar(:,2),'o');
title("Ilustrarea dependenței y^⋆(F^⋆_m) la momentul final t_f")
legend("theta1","theta2");
xlabel("F⋆ (N)");
ylabel("Theta ⋆(rad)")

%% Cerinta 8 
p1=polyfit(Fstar,ystar(:,1),8)
p2=polyfit(Fstar,ystar(:,2),8)

Fstar_grid=Fstar(1):0.01:Fstar(end);
y_interpolat1=polyval(p1,Fstar_grid);
y_interpolat2=polyval(p2,Fstar_grid);

plot(Fstar,ystar(:,1),'x',Fstar,ystar(:,2),'o');
hold on;
plot(Fstar_grid,y_interpolat1,'b--');
plot(Fstar_grid,y_interpolat2,'r--');
legend("theta1","theta2","y_interpolat1","y_interpolat2")
title("Gasirea polinomului de aproximare")
xlabel("F⋆ (N)")
ylabel("Theta⋆ (rad)")
hold off;


%% Cerinta 9 

alpha=randn(1,100)/10;

figure;
hold on;
for i=1:length(alpha)
    x0 = [theta1_init; (1+alpha(i))*theta2_init; 0; 0];
    out=sim("model_simulink");
    plot(out.tout,out.theta1.Data,'b',out.tout,out.theta2.Data,'r')
end
title("Incertitudine multiplicativă")
xlabel("timp (s)")
ylabel("theta (rad)")
legend("theta1","theta2");
hold off


%% Cerinta 10 

alpha=randn(1,100)*5;

figure;
hold on;
for i=1:length(alpha)
    x0 = [alpha(i)+theta1_init; theta2_init; 0; 0];
    out=sim("model_simulink");
    plot(out.tout,out.theta1.Data,'b',out.tout,out.theta2.Data,'r')
end
title("Incertitudine aditivă")
xlabel("timp (s)")
ylabel("theta (rad)")
legend("theta1","theta2");
hold off



%% Cerinta 11 

u=randn(1,length(t));
F=timeseries(u,t);
out=sim("model_simulink");

figure
set(gcf,"Position",[100 100 1200 600])
subplot(1,2,1)
dd=diff(diff(out.theta1.Data));
plot(out.tout(1:end-2),dd)
xlabel("timp(s)")
ylabel("ddtheta (rad/s^2)")
hold on;
yline(mean(dd),"r")
yline(3*std(dd),"g")
yline(-3*std(dd),"g")
legend("ddtheta","media","deviatia*3")
ylim([-1,1])
title("A doua derivata discreta a ieșirii")

hold off

subplot(1,2,2)
plot(t,u)
xlabel("timp (s)")
ylabel("F(N)")
hold on;
yline(mean(u),"r")
yline(3*std(u),"g")
yline(-3*std(u),"g")
ylim([-4,4])
title("Semnalul exogen F")
legend("F","media","deviatia*3")
hold off



