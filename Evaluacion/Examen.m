clear;
close all;

% Datos de entrada
v = [1.432 0.0 1.432 0.0 1.432 0.0 1.432 0.0 1.432 0.0, ]; 
w = [0.0 2.513 0.0 2.513 0.0 2.513 0.0 2.513 0.0 2.513, ]; 
dt = 1.0;  % Tiempo por paso

% Estado inicial
x_sim(1) = 0;
y_sim(1) = 0;
theta_sim(1) = 0;

% Simulación del movimiento
for i = 1:length(v)
    if w(i) == 0
        % Movimiento recto
        x_sim(i+1) = x_sim(i) + v(i) * dt * cos(theta_sim(i));
        y_sim(i+1) = y_sim(i) + v(i) * dt * sin(theta_sim(i));
        theta_sim(i+1) = theta_sim(i);
    else
        % Giro en el lugar
        x_sim(i+1) = x_sim(i);
        y_sim(i+1) = y_sim(i);
        theta_sim(i+1) = theta_sim(i) + w(i) * dt;
    end
    i
    fprintf('x = %.3f m\ny = %.3f m\ntheta = %.3f rad/s\n', x_sim(i+1), y_sim(i+1), theta_sim(i+1));
end



% Animación
x1 = x_sim;
y1 = y_sim;
phi = theta_sim;
hx = x_sim;
hy = y_sim;
ts = 0.1;         % Tiempo entre cuadros
N = length(x1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VISUALIZACIÓN Y GRÁFICAS 3D %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Configuración de escena
scene = figure;
set(scene,'Color','white');
set(gca,'FontWeight','bold');
sizeScreen = get(0,'ScreenSize');
set(scene,'position',sizeScreen);
camlight('headlight');
axis equal;
grid on;
box on;
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
view([15 15]);
axis([-6 8 -12 8 0 2]);  % Ajustado al dominio y rango de la función

% Inicialización del modelo del robot
scale = 4;
MobileRobot_5;  
H1 = MobilePlot_4(x1(1), y1(1), phi(1), scale); hold on;
H2 = plot3(hx(1), hy(1), 0, 'r', 'lineWidth', 2);

% Animación
step = 1;
for k = 1:step:N
    delete(H1);    
    delete(H2);

    H1 = MobilePlot_4(x1(k), y1(k), phi(k), scale);
    H2 = plot3(hx(1:k), hy(1:k), zeros(1,k), 'r', 'lineWidth', 2);

    pause(ts);
end
