clc; clear; close all;

% Trayectoria 1
%x = linspace(0, 5, 1000);            % Espacio de tiempo o paso en x
%y = 2 * sin(x.^2);                  % Trayectoria deseada
%dt = x(2) - x(1);                   % Paso de tiempo equivalente

%Trayectoria 2
path = linspace(0,2*pi,1000);
r = 4;
x = r*cos(path);
y = r*sin(path);
dt = path(2) - path(1)

% Funcion 3
%x = linspace(-6, 6, 1000);  % Rango de x
%y = zeros(size(x));        % Inicializamos y

%for i = 1:length(x)
%    if x(i) <= -1
%        y(i) = 2 * x(i);
%    elseif x(i) < 1
%        y(i) = 2 * x(i) + 1;
%    elseif x(i) < 4
%        y(i) = -x(i) + 4;
%    else
%        y(i) = x(i) - 1;
%    end
%end

%dt = x(2) - x(1);

dx = gradient(x, dt);
dy = gradient(y, dt);

v = sqrt(dx.^2 + dy.^2);         % Velocidad lineal
theta = atan2(dy, dx);           % Orientación
omega = gradient(theta, dt);     % Velocidad angular


x_sim = zeros(1, length(x));
y_sim = zeros(1, length(x));
theta_sim = zeros(1, length(x));

x_sim(1) = x(1);
y_sim(1) = y(1);
theta_sim(1) = theta(1);


for i = 2:length(x)
    x_sim(i) = x_sim(i-1) + v(i-1) * cos(theta_sim(i-1)) * dt;
    y_sim(i) = y_sim(i-1) + v(i-1) * sin(theta_sim(i-1)) * dt;
    theta_sim(i) = theta_sim(i-1) + omega(i-1) * dt;
end

% Animacion
x1 = x_sim;
y1 = y_sim;
phi = theta_sim;
hx = x_sim;
hy = y_sim;
ts = 0.01;         % Tiempo entre cuadros
N = length(x1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VISUALIZACIÓN Y GRÁFICAS 3D %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Configuración de escena
scene=figure;
set(scene,'Color','white');
set(gca,'FontWeight','bold');
sizeScreen=get(0,'ScreenSize');
set(scene,'position',sizeScreen);
camlight('headlight');
axis equal;
grid on;
box on;
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
view([15 15]);
axis([-6 8 -12 8 0 2]);  % Ajustado al dominio y rango de la función


scale = 4;
MobileRobot_5;  % Asume que carga el modelo o configuración del robot
H1 = MobilePlot_4(x1(1), y1(1), phi(1), scale); hold on;


H2 = plot3(hx(1), hy(1), 0, 'r', 'lineWidth', 2);


step = 5;  % Ajusta para velocidad de animación

for k = 1:step:N
    delete(H1);    
    delete(H2);

    H1 = MobilePlot_4(x1(k), y1(k), phi(k), scale);
    H2 = plot3(hx(1:k), hy(1:k), zeros(1,k), 'r', 'lineWidth', 2);

    pause(ts);
end