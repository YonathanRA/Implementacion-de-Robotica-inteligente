% Create waypoints
op = 3;

if  op == 1
    waypoints = [2,8; 
                 8, 8;
                 8, 2;
                 2, 2];
elseif op ==2
    waypoints = [1 2;
                 8 9;
                 6 7;
                 3 9;
                 1 7;
                 6 2;
                 8 3];
elseif op ==3
    waypoints = [1 1; 
                 1 4;
                 1 7;
                 1 10;
                 4 4;
                 4 7;
                 4 10;
                 7 1;
                 7 4;
                 7 7;
                 10 1;
                 10 4;
                 10 7;
                 10 10];
elseif op == 4
    waypoints = [1 2; 
                 2 10;
                 11 8;
                 8 2;
                 8, 8;
                 1 2];
end

%% Simulation setup
% Define Vehicle
R = 0.1;                        % Wheel radius [m]
L = 0.5;                        % Wheelbase [m]
dd = DifferentialDrive(R,L);

% Sample time and time array
sampleTime = 0.05;              % Sample time [s]
tVec = 0:sampleTime:190;        % Time array

% Initial conditions
initPose = [waypoints(1,1);waypoints(1,2); 0];            % Initial pose (x y theta)
pose = zeros(3,numel(tVec));   % Pose matrix
pose(:,1) = initPose;


% Load map

%complexMap       41x52                2132  logical              
%emptyMap         26x27                 702  logical              
%simpleMap        26x27                 702  logical              
%ternaryMap      501x501            2008008  double  

close all
load exampleMap
%load complexMap

% Create lidar sensor
lidar = LidarSensor;
lidar.sensorOffset = [0,0];
lidar.scanAngles = linspace(-pi,pi,90);%51
lidar.maxRange = 5;%5


% Create visualizer
viz = Visualizer2D;
viz.hasWaypoints = true;
viz.mapName = 'map';
attachLidarSensor(viz,lidar);

%% Path planning and following

% Pure Pursuit Controller
controller = controllerPurePursuit;
controller.Waypoints = waypoints;
controller.LookaheadDistance = 1;%0.5
controller.DesiredLinearVelocity = 0.45; %0.75
controller.MaxAngularVelocity = 8;

% Vector Field Histogram (VFH) for obstacle avoidance
vfh = controllerVFH;
vfh.DistanceLimits = [0.08 5]; %0.05 3
vfh.NumAngularSectors = 38; %36
vfh.HistogramThresholds = [10 15]; % 5y 10
vfh.RobotRadius = L*6/8;
vfh.SafetyDistance = L/2;
vfh.MinTurningRadius = 0.65;%0.25

%% Simulation loop
r = rateControl(1/sampleTime);
for idx = 2:numel(tVec) 
    
    % Get the sensor readings
    curPose = pose(:,idx-1);
    ranges = lidar(curPose);
        
    % Run the path following and obstacle avoidance algorithms
    [vRef,wRef,lookAheadPt] = controller(curPose);
    targetDir = atan2(lookAheadPt(2)-curPose(2),lookAheadPt(1)-curPose(1)) - curPose(3);
    steerDir = vfh(ranges,lidar.scanAngles,targetDir);    
    if ~isnan(steerDir) && abs(steerDir-targetDir) > 0.1
        wRef = 0.5*steerDir;
    end
    
    % Control the robot
    velB = [vRef;0;wRef];                   % Body velocities [vx;vy;w]
    vel = bodyToWorld(velB,curPose);  % Convert from body to world
    
    % Perform forward discrete integration step
    pose(:,idx) = curPose + vel*sampleTime; 
    
    % Update visualization
    viz(pose(:,idx),waypoints,ranges)
    waitfor(r);
end