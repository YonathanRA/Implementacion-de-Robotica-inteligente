clear 
close all
clc

% Control de tipo de control
% Lazo abierto 1
% Lazo Cerrado 2
% COntrol Velocidad 3

op = 3;

% Vectores de coordenadas
x_vec = [ -9 , -7 , -1  , 3 ,  5 , -9  ];
y_vec = [  10 , 1 , -10 , 1 , 10 ,  10 ];

if op ==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Lazo Abierto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parámetros del movimiento
ts = 0.01;       % Tiempo de muestreo
v = 10;         % Velocidad lineal
w = 5;         % Velocidad angular

% Inicialización del robot
x = x_vec(1); y = y_vec(1); theta = 0; % Posición y orientación inicial
x1 = []; y1 = []; phi = [];            % Trayectoria simulada
hx = []; hy = [];                      % Historia del centro del robot

for i = 1:(length(x_vec)-1)
    x_goal = x_vec(i+1);
    y_goal = y_vec(i+1);
    
    % ROTACIÓN EN EL LUGAR
    % Calcular el ángulo al objetivo
    angle_to_goal = atan2(y_goal - y, x_goal - x);
    
    % Calcular la diferencia de ángulos más corta
    angle_diff = angle_to_goal - theta;
    
    % Normalizar la diferencia de ángulos para estar en el rango [-pi, pi]
    angle_diff = mod(angle_diff + pi, 2*pi) - pi;
    
    % Calcular el tiempo necesario para rotar
    t_rot = abs(angle_diff) / w;  % Tiempo en segundos para llegar al ángulo
    
    % Actualizar la orientación del robot
    theta = theta + angle_diff;
    
    % Registrar la trayectoria
    x1(end+1) = x;
    y1(end+1) = y;
    phi(end+1) = theta;
    hx(end+1) = x;
    hy(end+1) = y;

    % Calcular la distancia al objetivo
    dist = sqrt((x_goal - x)^2 + (y_goal - y)^2);
    
    % Calcular el tiempo necesario para moverse hasta el objetivo
    t_trans = dist / v;  % Tiempo en segundos para cubrir la distancia
    
   % Movimiento lineal paso a paso
    n_steps = ceil(t_trans / ts); % Número de pasos
    for j = 1:n_steps
        x = x + v * cos(theta) * ts;
        y = y + v * sin(theta) * ts;
    
        % Registrar trayectoria paso a paso
        x1(end+1) = x;
        y1(end+1) = y;
        phi(end+1) = theta;
        hx(end+1) = x;
        hy(end+1) = y;
    end
end

% Tiempo total de simulación
k_global = length(x1);

elseif op == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Lazo Cerrado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tiempo de simulación y muestreo
tf = 10*length(x_vec);          % Tiempo total (s)
ts = 0.1;                       % Tiempo de muestreo (s)
t = 0:ts:tf;
N = length(t);

% Condiciones iniciales
x1(1) = x_vec(1);
y1(1) = y_vec(1);
phi(1) = 0;

% Punto de control
hx(1) = x1(1);
hy(1) = y1(1);

% Inicializamos vectores de velocidad
v = zeros(1, N);
w = zeros(1, N);
Error = zeros(1, N);

k_global = 1; % índice global para recorrer todo el vector de tiempo

% Bucle principal para recorrer todos los puntos de x_vec e y_vec
for i = 1:length(x_vec)

    hxd = x_vec(i);
    hyd = y_vec(i);

    while k_global < N
        % Errores de posición
        hxe(k_global) = hxd - hx(k_global);
        hye(k_global) = hyd - hy(k_global);
        Error(k_global) = sqrt(hxe(k_global)^2 + hye(k_global)^2);

        % Condición de llegada al punto
        if Error(k_global) < 0.17
            break;
        end

        % Jacobiano
        J = [cos(phi(k_global)) -sin(phi(k_global));
             sin(phi(k_global))  cos(phi(k_global))];

        % Ganancias
        K = [1.5 0;
             0 1.5];

        % Ley de control
        he = [hxe(k_global); hye(k_global)];
        qpRef = pinv(J) * K * he;
        v(k_global) = qpRef(1);
        w(k_global) = qpRef(2);

        % Actualización del estado
        phi(k_global + 1) = phi(k_global) + w(k_global) * ts;
        xp1 = v(k_global) * cos(phi(k_global));
        yp1 = v(k_global) * sin(phi(k_global));
        x1(k_global + 1) = x1(k_global) + xp1 * ts;
        y1(k_global + 1) = y1(k_global) + yp1 * ts;
        hx(k_global + 1) = x1(k_global + 1);
        hy(k_global + 1) = y1(k_global + 1);

        k_global = k_global + 1;
    end

end

end

if op == 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Control de Velocidades
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 TIEMPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf=20;             % Tiempo de simulación en segundos (s)
ts=0.05;            % Tiempo de muestreo en segundos (s)
t=0:ts:tf;         % Vector de tiempo
N= length(t);      % Muestras


%2 CONDICIONES INICIALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damos valores a nuestro punto inicial de posición y orientación
x1(1)= -9;  %Posición inicial eje x
y1(1)= 10;  %Posición inicial eje y
phi(1)= -pi/2; %Orientación inicial del robot 

%Igualamos el punto de control con las proyecciones X1 y Y1 por su
%coincidencia
hx(1)= x1(1);       % Posición del punto de control en el eje (X) metros (m)
hy(1)= y1(1);       % Posición del punto de control en el eje (Y) metros (m)

%3 TRAYECTORIA DESEADA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_lado = tf/5;

% Vectores de salida
hxd = zeros(1,N);
hyd = zeros(1,N);
hxdp = zeros(1,N);
hydp = zeros(1,N);

% Para cada tramo
for i = 1:5
    idx_start = floor((i-1)*T_lado/ts) + 1;
    idx_end = min(floor(i*T_lado/ts), N);
    idx = idx_start:idx_end;
    len = length(idx);

    x0 = x_vec(i);     x1_ = x_vec(i+1);
    y0 = y_vec(i);     y1_ = y_vec(i+1);

    hxd(idx) = linspace(x0, x1_, len);
    hyd(idx) = linspace(y0, y1_, len);

    vx_lado = (x1_ - x0) / T_lado;
    vy_lado = (y1_ - y0) / T_lado;
    hxdp(idx) = vx_lado * ones(1, len);
    hydp(idx) = vy_lado * ones(1, len);
end

%4 CONTROL, BUCLE DE SIMULACION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N 

    %a)Errores de control (Aqui la posición deseada ya no es constante,
    % varia con el tiempo)
    hxe(k)=hxd(k)-hx(k);
    hye(k)=hyd(k)-hy(k);
    
    %Matriz de error
    he= [hxe(k);hye(k)];
    %Magnitud del error de posición
    Error(k)= sqrt(hxe(k)^2 +hye(k)^2);

    %b)Matriz Jacobiana
    J=[cos(phi(k)) -sin(phi(k));... %Matriz de rotación en 2D
       sin(phi(k)) cos(phi(k))];

    %c)Matriz de Ganancias
    K=[5 0;...
       0 5];
    
    %d)Velocidades deseadas
    hdp=[hxdp(k);hydp(k)];

    %e)Ley de Control:Agregamos las velocidades deseadas
    qpRef= pinv(J)*(hdp + K*he);

    v(k)= qpRef(1);   %Velocidad lineal de entrada al robot 
    w(k)= qpRef(2);   %Velocidad angular de entrada al robot 


%5 APLICACIÓN DE CONTROL AL ROBOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Aplico la integral a la velocidad angular para obtener el angulo "phi" de la orientación
    phi(k+1)=phi(k)+w(k)*ts; % Integral numérica (método de Euler)
           
   %%%%%%%%%%%%%%%%%%%%% MODELO CINEMATICO %%%%%%%%%%%%%%%%%%%%%%%%%
    
    xp1=v(k)*cos(phi(k)); 
    yp1=v(k)*sin(phi(k));
 
    %Aplico la integral a la velocidad lineal para obtener las cordenadas
    %"x1" y "y1" de la posición
    x1(k+1)=x1(k)+ ts*xp1; % Integral numérica (método de Euler)
    y1(k+1)=y1(k)+ ts*yp1; % Integral numérica (método de Euler)

    % Posicion del robot con respecto al punto de control
    hx(k+1)=x1(k+1); 
    hy(k+1)=y1(k+1);
     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULACION VIRTUAL 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a) Configuracion de escena

scene=figure;  % Crear figura (Escena)
set(scene,'Color','white'); % Color del fondo de la escena
set(gca,'FontWeight','bold') ;% Negrilla en los ejes y etiquetas
sizeScreen=get(0,'ScreenSize'); % Retorna el tamaño de la pantalla del computador
set(scene,'position',sizeScreen); % Configurar tamaño de la figura
camlight('headlight'); % Luz para la escena
axis equal; % Establece la relación de aspecto para que las unidades de datos sean las mismas en todas las direcciones.
grid on; % Mostrar líneas de cuadrícula en los ejes
box on; % Mostrar contorno de ejes
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); % Etiqueta de los eje

view([-0.1 90]); % Orientacion de la figura
axis([-15 15 -15 15 0 1]); % Ingresar limites minimos y maximos en los ejes x y z [minX maxX minY maxY minZ maxZ]

scale = 4;
MobileRobot_5;
H1=MobilePlot_4(x1(1),y1(1),phi(1),scale);hold on;

% c) Graficar Trayectorias
H2=plot3(hx(1),hy(1),0,'r','lineWidth',2);
H3=plot3(hxd,hyd,zeros(1,N),'g','lineWidth',2);

step=1; % pasos para simulacion

for k=1:step:N

    delete(H1);    
    delete(H2);
    
    H1=MobilePlot_4(x1(k),y1(k),phi(k),scale);
    H2=plot3(hx(1:k),hy(1:k),zeros(1,k),'r','lineWidth',2);
    
    pause(ts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graficas %%%%%%%%%%%%%%%%%%%%%%%%%%%%
graph=figure;  % Crear figura (Escena)
set(graph,'position',sizeScreen); % Congigurar tamaño de la figura
subplot(311)
plot(t,v,'b','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('m/s'),legend('Velocidad Lineal (v)');
subplot(312)
plot(t,w,'g','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[rad/s]'),legend('Velocidad Angular (w)');
subplot(313)
plot(t,Error,'r','LineWidth',2),grid('on'),xlabel('Tiempo [s]'),ylabel('[metros]'),legend('Error de posición (m)');


else
% Debido a que Lazo cerrado y Lazo abierto usan la misma forma de graficar
% se guarda en un else, para asi poder manejarlo distinto a Control de
% Velocidad
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULACIÓN VIRTUAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configuración de escena
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
view([0 90]);

% Límites automáticos
xmin_lim = min(x_vec) - 2;
xmax_lim = max(x_vec) + 2;
ymin_lim = min(y_vec) - 2;
ymax_lim = max(y_vec) + 2;
axis([xmin_lim xmax_lim ymin_lim ymax_lim 0 1]);

% Gráfica inicial
scale = 4;
MobileRobot_5; 
H1 = MobilePlot_4(x1(1), y1(1), phi(1), scale); hold on;
H2 = plot3(hx(1), hy(1), 0, 'r', 'lineWidth', 2);
H3 = plot3(x_vec, y_vec, 0 * ones(size(x_vec)), 'bo', 'lineWidth', 2);
H4 = plot3(hx(1), hy(1), 0, 'go', 'lineWidth', 2);

% Animación
step = 1;
for k = 1:step:k_global
    delete(H1);
    delete(H2);
    H1 = MobilePlot_4(x1(k), y1(k), phi(k), scale);
    H2 = plot3(hx(1:k), hy(1:k), zeros(1,k), 'r', 'lineWidth', 2);
    pause(ts);
end

end