clear; close all;

% Parámetros del robot
R = 0.05;    % Radio de las ruedas (m)
L = 0.18;    % Distancia entre ruedas (m)
Rc = 15;     % Radio del círculo deseado (m)
v_const = 0.65;        % Velocidad lineal constante (m/s)
w_const = v_const / Rc;  % Velocidad angular para un círculo de radio Rc
dt = 1;      % Intervalo de tiempo (s)
T_total =  2 * pi * Rc / v_const;  % Tiempo total para una vuelta completa
n = round(T_total / dt);          % Número de pasos

% Inicialización
v = v_const * ones(n,1);
w = w_const * ones(n,1);
x_sim = zeros(n+1,1);
y_sim = zeros(n+1,1);
theta_sim = zeros(n+1,1);

% Centro en (0,0), iniciar en (0, -Rc)
x_sim(1) = 0;
y_sim(1) = -Rc;
theta_sim(1) = 0;

% Simulación de trayectoria
for i = 1:n
    x_sim(i+1) = x_sim(i) + v(i) * dt * cos(theta_sim(i));
    y_sim(i+1) = y_sim(i) + v(i) * dt * sin(theta_sim(i));
    theta_sim(i+1) = theta_sim(i) + w(i) * dt;
end

% Guardar CSV
vel_data = table((1:n)', v, w, 'VariableNames', {'Paso', 'v_m_s', 'w_rad_s'});
writetable(vel_data, 'velocidades_robot.csv');

% Animación
x1 = x_sim;
y1 = y_sim;
phi = theta_sim;
hx = x_sim;
hy = y_sim;
ts = 0.01;
N = length(x1);

% Escena 3D
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
xlim([-Rc-2 Rc+2]);
ylim([-Rc-2 Rc+2]);
zlim([0 2]);

% Inicializar robot
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
fprintf('x = %.3f m\ny = %.3f m\n', x1(end), y1(end));

text(x1(end), y1(end), 0.5, sprintf('x=%.2f m\\ny=%.2f m', x1(end), y1(end)), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k');
