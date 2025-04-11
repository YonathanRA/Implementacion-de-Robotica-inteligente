% Parámetros

% Trayectoria 1
%x = linspace(0, 5, 1000);            % Espacio de tiempo o paso en x
%y = 2 * sin(x.^2);                  % Trayectoria deseada
%dt = x(2) - x(1);                   % Paso de tiempo equivalente

%Trayectoria 2
%path = linspace(0,2*pi,1000);
%r = 4;
%x = r*cos(path);
%y = r*sin(path);
%dt = path(2) - path(1)

%Trayectoria 3
x = linspace(-6, 6, 1000000); 
y = zeros(size(x));

for i = 1:length(x)
    if x(i) <= -1
        y(i) = 2 * x(i);
    elseif x(i) < 1
        y(i) = 2 * x(i) + 1;
    elseif x(i) < 4
        y(i) = -x(i) + 4;
    else
        y(i) = x(i) - 1;
    end
end

dt = x(2)-x(1);

% Derivadas para obtener velocidad y orientación
dx = gradient(x, dt);
dy = gradient(y, dt);

% Magnitud de la velocidad (lineal)
v = sqrt(dx.^2 + dy.^2);

% Ángulo de orientación (theta)
theta = atan2(dy, dx);

% Derivada de theta para obtener velocidad angular omega
omega = gradient(theta, dt);

% Simular trayectoria con velocidades (la integración)
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

% Plot para comparar
figure;
plot(x, y, 'b--', 'LineWidth', 2); hold on;
plot(x_sim, y_sim, 'r', 'LineWidth', 2);
legend('Trayectoria deseada', 'Trayectoria del robot');
xlabel('x'); ylabel('y');
title('Seguimiento de trayectoria en lazo abierto');
grid on;