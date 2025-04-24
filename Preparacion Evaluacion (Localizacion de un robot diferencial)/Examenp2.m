% Leer tabla
entrada = readtable('entrada.csv');

% Parámetros del robot
R = 0.1;    % Radio de las ruedas (m)
L = 0.4;    % Distancia entre ruedas (m)
dt = 1;     % intervalo de tiempo (s)

n = height(entrada);
v = zeros(n,1);
w = zeros(n,1);
x = zeros(n+1,1);
y = zeros(n+1,1);
theta = zeros(n+1,1); 

% Cálculo de v, w y trayectoria
for i = 1:n
    omega_R = entrada.omega_R(i);
    omega_L = entrada.omega_L(i);
    
    v(i) = (R/2) * (omega_R + omega_L);
    w(i) = (R/L) * (omega_R - omega_L);
    
    x(i+1) = x(i) + v(i) * dt * cos(theta(i));
    y(i+1) = y(i) + v(i) * dt * sin(theta(i));
    theta(i+1) = theta(i) + w(i) * dt;
end

% Crear tabla
salida = table((0:n)', [v; NaN], [w; NaN], x, y, rad2deg(theta), ...
    'VariableNames', {'t_s', 'v_mps', 'w_radps', 'x_m', 'y_m', 'theta_deg'});

% Guardar archivo CSV
writetable(salida, 'salida.csv');

% Graficar
figure;
plot(x, y, 'b-o', 'LineWidth', 2); hold on;
plot(x(1), y(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x(end), y(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on; axis equal;
legend('Trayectoria', 'Inicio', 'Fin');