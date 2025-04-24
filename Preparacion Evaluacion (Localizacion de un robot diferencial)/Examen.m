% Datos de entrada
v = [1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0];  % velocidad lineal
w = [0.0 pi/3 0.0 pi/3 0.0 pi/3 0.0 pi/3 0.0 pi/3 0.0 pi/3];  % velocidad angular
dt = 1.0;  % tiempo por paso

% Estado inicial
x(1) = -1;
y(1) = -5;
theta(1) = 0;

% Simulación
for i = 1:length(v)
    if w(i) == 0
        % Movimiento recto
        x(i+1) = x(i) + v(i) * dt * cos(theta(i));
        y(i+1) = y(i) + v(i) * dt * sin(theta(i));
        theta(i+1) = theta(i);
    else
        % Giro en el lugar
        x(i+1) = x(i);
        y(i+1) = y(i);
        theta(i+1) = theta(i) + w(i) * dt;
    end
end

% Graficar la trayectoria
figure;
plot(x, y, 'o-', 'LineWidth', 2); hold on;
grid on;
xlabel('x (m)');
ylabel('y (m)');
title('Trayectoria del robot diferencial');
axis equal;

% Marcar punto de inicio y fin
plot(x(1), y(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');  % inicio
plot(x(end), y(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');  % fin

legend('Trayectoria', 'Inicio', 'Fin');
