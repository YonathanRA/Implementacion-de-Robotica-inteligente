clear;
close all;

% Leer tabla
entrada = readtable('entrada.csv');

% Parámetros del robot
R = 0.05;    % Radio de las ruedas (m)
L = 0.18;    % Distancia entre ruedas (m)
dt = 1;      % Intervalo de tiempo (s)

n = height(entrada);
v = zeros(n,1);
w = zeros(n,1);
x_sim = zeros(n+1,1);
y_sim = zeros(n+1,1);
theta_sim = zeros(n+1,1); 

% Cálculo de v, w y trayectoria
for i = 1:n
    omega_R = entrada.omega_R(i);
    omega_L = entrada.omega_L(i);
    
    v(i) = (R/2) * (omega_R + omega_L);
    w(i) = (R/L) * (omega_R - omega_L);
    
    x_sim(i+1) = x_sim(i) + v(i) * dt * cos(theta_sim(i));
    y_sim(i+1) = y_sim(i) + v(i) * dt * sin(theta_sim(i));
    theta_sim(i+1) = theta_sim(i) + w(i) * dt;
end

% Crear tabla
salida = table((0:n)', [v; NaN], [w; NaN], x_sim, y_sim, rad2deg(theta_sim), ...
    'VariableNames', {'t_s', 'v_mps', 'w_radps', 'x_m', 'y_m', 'theta_deg'});

% Guardar archivo CSV
writetable(salida, 'salida.csv');


% Animación
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

% Límites automáticos basados en trayectoria
xlim([min(x1)-1 max(x1)+1])
ylim([min(y1)-1 max(y1)+1])
zlim([0 2])

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

% Mostrar posición final
fprintf('Posición final del robot:\n');
fprintf('x = %.3f m | y = %.3f m', x1(end), y1(end));

text(x1(end), y1(end), 0.5, sprintf('x=%.2f m\\ny=%.2f m', x1(end), y1(end)), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k');
