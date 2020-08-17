%% Thesis Sonar Model  
%% Andi Rascon, Woen-Sug Choi


clear;clc;close all;

%% Input parameters 

sourcePos = [0;5;8];
sourceNormal = [0.7065;0.48;-0.52];
meshResolution = 30;
nRays = 5;  % # of rays (ONLY ODD NUMBER working on angle calculation )
nBeams = 50; % # number of Beams

%-------------------------------------
% Given Input parameters
%-------------------------------------
% Sonar properties
% BlueView P900-90 Modified
sonarFreq = 900E3; % Sonar frquency
maxRange = 60; % Maximum unambiguous range
rangeRes = 0.025; % Range resolution (1 inch = 0.0254 meter)
beamWidth = 30/360*2*pi(); % [rad]
fieldOfView = 45; % [degree]
bandwidth = 2.95e4; % [Hz]
S_f = exp(-(f-sonarFreq).^2*bandwidth^2*pi^2); %transmit spectrum, frequency domain
%-------------------------------------
% Computational input parameters
%-------------------------------------
% Sonar transducer eletrical properties
SL = 220; % [dB re 1uPa at 1m]

% Environmental properties
TSvalue = 25; % further:academic.oup.com/icesjms/article/71/4/882/666881
soundSpeed = 1500.0; % [m/s]


%-------------------------------------
% Bathymetry Grids
%-------------------------------------

% create cylinder
r = 5;  % radius
h = 20; % height
resolution = 10; % resolution 
% Make triangle geometry
[X,Y,Z] = cylinder(r,resolution);
Z = Z*h;
% plot
surf(X,Y,Z);

% Surface Reflectivity
mu = 10e-4;


[xMesh,yMesh] = meshgrid(1:meshResolution,1:meshResolution);
zMesh = peaks(meshResolution); 
% zMesh = zeros(meshResolution);
gridConnectivity = delaunay(xMesh,yMesh);
nTargets = length(gridConnectivity);
% xMesh(1:21,1:2) = X.';
% yMesh(1:21,1:2) = Y.';
% zMesh(1:21,1:2) = Z.';
targetPoints = [xMesh(:), yMesh(:), zMesh(:)];
% to plot :
% trimesh(gridConnectivity,X,Y,Z);
% mesh(X,Y,Z)
% hold on;

% Rearrange array for gpu ray tracing
P0 = zeros(nTargets,3);
P1 = zeros(nTargets,3);
P2 = zeros(nTargets,3);
for i = 1:nTargets
    P0(i,:) = targetPoints(gridConnectivity(i,1),:);
    P1(i,:) = targetPoints(gridConnectivity(i,2),:);
    P2(i,:) = targetPoints(gridConnectivity(i,3),:);
end


%% Transmitted acoustics

% Transmitted Waveform properties
prf = soundSpeed/(2*maxRange);            % Pulse repetition frequency
pulse_width = 2*rangeRes/soundSpeed;      % Pulse width
pulse_bw = 1/pulse_width;                 % Pulse bandwidth
fs = 2*pulse_bw;                          % Sampling rate

% Generate rectangular pulse
Pulse = zeros(fs/prf,1); 
for i=1:fs*pulse_width; 
    Pulse(i) = 1; 
end


%% Beam pattern function
Lrambda = 0.884/beamWidth;
BeamPattern = @(theta) abs(sin(pi()*Lrambda.*sin(theta))./(pi()*Lrambda.*sin(theta))).^2;
% To plot :
% theta = -beamWidth:0.01:beamWidth;
% beamPatternPlot = BeamPattern(theta);
% polarplot(theta,beamPatternPlot); thetalim([-beamWidth/2/pi()*360 beamWidth/2/pi()*360])


%% Sonar simulation

timeDomain = (0:size(Pulse,1)-1)/fs;
SignalBeam = zeros(size(Pulse,1),nBeams);

% Calculate beam angles
centerAzimuthAngle = atan(sourceNormal(2)/sourceNormal(1));
leftAzimuthAngle = centerAzimuthAngle-(fieldOfView/360*2*pi())/2;
rightAzimuthAngle = centerAzimuthAngle+(fieldOfView/360*2*pi())/2;
beamAzimuthAngles = linspace(leftAzimuthAngle,rightAzimuthAngle,nBeams);
centerElevationAngle = atan(sqrt(sourceNormal(1)^2+sourceNormal(2)^2)/sourceNormal(3));
topElevationAngle = centerElevationAngle+beamWidth/2;
bottomElevationAngle = centerElevationAngle-beamWidth/2;
elevationAngles = linspace(bottomElevationAngle,topElevationAngle,nRays);

% Subplot for rays intersecting grid
subplot(1,2,1)
trimesh(gridConnectivity,xMesh,yMesh,zMesh); hold on; 

% create frequency vector
fmin = sonarFreq - bandwidth/2;
fmax = sonarFreq + bandwidth/2;
f = linspace(fmin,fmax,1e5);


% Loop for Beams
f = waitbar(0,'Calculating Beams');
for i = 1:nBeams
    % --- Loop for rays to calculate ray-tracing data --- %
    % Calculate center, top, bottom degree of a beam
    azimuthAngle = beamAzimuthAngles(i);
    
    % Ray loop
    rayData = zeros(5,nRays); % rayData array for distance,grid number,normalVector
    for k = 1:nRays
        % Calculate ray directions in a beam
        rayDirections = [sin(elevationAngles(k))*cos(azimuthAngle);sin(elevationAngles(k))*sin(azimuthAngle);cos(elevationAngles(k))];
        % --- Use gpu ray tracing external function --- %
        % mathworks.com/matlabcentral/fileexchange/49670-hardware-accelerated-ray-triangle-intersection
        % Find out which grid intersect and its distance
        or = ones(nTargets,3);
        or(:,1) = sourcePos(1); or(:,2) = sourcePos(2); or(:,3) = sourcePos(3);
        D = ones(nTargets,3);
        D(:,1) = rayDirections(1); D(:,2) = rayDirections(2); D(:,3) = rayDirections(3);
        gD = gpuArray(D);
        [dist, flag] = arrayfun(@rayTriGPU,P0(:,1)',P0(:,2)',P0(:,3)', ...
            P1(:,1)', P1(:,2)', P1(:,3)', ...
            P2(:,1)', P2(:,2)', P2(:,3)', ...
            or(:,1), or(:,2), or(:,3), ...
            gD(:,1),gD(:,2),gD(:,3));
        distances = gather(dist);
        % Find out index of intersecting grid
        for j = 1:nTargets
            if isnan(distances(1,j)) == false
                rayData(1,k) = abs(distances(1,j)); % save distance
                rayData(2,k) = j; % save grid number
                % calculate normal
                vector1 = [P1(j,1)-P0(j,1) P1(j,2)-P0(j,2) P1(j,3)-P0(j,3)];
                vector2 = [P2(j,1)-P0(j,1) P2(j,2)-P0(j,2) P2(j,3)-P0(j,3)];
                normalVector = cross(vector1,vector2);
                normalVector = normalVector/norm(normalVector);
                rayData(3,k) = normalVector(1);
                rayData(4,k) = normalVector(2);
                rayData(5,k) = normalVector(3);
            end
        end
    end
    
    
    % --- Acoustics calculation --- %
    Signal = zeros(size(Pulse,1),nRays);
    
    % Loop over intersected targets
    rayEL = 0;
    amplitude = 0;

    
    for k = 1:nRays
        % Loop for ray intersected targets to calculate sonar equations
        % EL = SL-2TL+TS
        % Ray tracing data allocation
        rayDirections = [sin(elevationAngles(k))*cos(azimuthAngle);sin(elevationAngles(k))*sin(azimuthAngle);cos(elevationAngles(k))];
        distance = rayData(1,k); gridNumber = rayData(2,k);
        gridNormal = [rayData(3,k); rayData(4,k); rayData(5,k)];
        % skip if ray tracing didnt intersect or went wrong
        if distance == 0
            continue;
        end
        
        %-- Source level --%
%         if elevationAngles(k)-centerElevationAngle == 0
%             raySL = SL;
%         else
%             raySL = SL - 10*log10(BeamPattern(elevationAngles(k)-centerElevationAngle));
%         end
        
        %-- Transmission path loss --% TL = 20log10(r) + alpha*r
%         alpha = @(freq) (3.0E-3 + (0.1*freq^2)/(1+freq^2)+(44*freq^2)/(4100+freq^2)+2.75E-4*freq^2);
%         absorptionCoeff = alpha(sonarFreq); % [db/km]
%         rayTL = 20*log10(distance); %;+ absorptionCoeff*distance/1000;
%         
%         %-- Target Strength --% simple cosine model
%         rayTS = TSvalue/cos(norm(gridNormal.*rayDirections)/(norm(gridNormal)*norm(rayDirections)));

x = rand;   % generate a random number, part of Gaussian noise  
y = rand;   % generate another random number, part of Gaussian noise
P_f(1) = 0;
Nscatters = 5;  % This should be the number of surfaces a ray intersects with. For now, using constant for simplicity

for j = 1:Nscatters
% scattering amplitude, includes overall Gaussian noise and target strength
% alpha needs to be the angle between incident ray and normal vector
amplitude(j) = ((x+i*y)/sqrt(2))*(sqrt(mu*cos(alpha)^2*distance^2*azimuthAngle*elevationAngle));

        
%         %-- Summation for final power level--%
%         rayEL = raySL - 2*rayTL + rayTS;
        
        %-- Signal Calculation --%
        rayRetTime = distance*2/soundSpeed;
        
  % Summation of Echo returned from a signal
        P_f(j) = P_f(j) + S_f*amplitdue*exp(i*2*pi*sonarFreq*rayRetTime/(2*distance))/(distance^2);
        
        
end
        Signal(:,k) = circshift(Pulse*rayEL,floor(rayRetTime*fs));
    end
    SignalBeam(:,i) = pulsint(Signal,'noncoherent');
    
    
    % Plot Rays and intersecting grids
    subplot(1,2,1)
    scatter3(sourcePos(1),sourcePos(2),sourcePos(3));hold on;
    for j = 1:nRays
        rayDirections = [sin(elevationAngles(j))*cos(azimuthAngle);sin(elevationAngles(j))*sin(azimuthAngle);cos(elevationAngles(j))];
        rayDraw = -rayData(1,j)*1.1*rayDirections+sourcePos;
        plot3([sourcePos(1) rayDraw(1)],[sourcePos(2) rayDraw(2)],[sourcePos(3) rayDraw(3)],'r');
        gridNumber = rayData(2,j);
        plot3([P0(gridNumber,1) P1(gridNumber,1) P2(gridNumber,1) P0(gridNumber,1)], ...
            [P0(gridNumber,2) P1(gridNumber,2) P2(gridNumber,2) P0(gridNumber,2)], ...
            [P0(gridNumber,3) P1(gridNumber,3) P2(gridNumber,3) P0(gridNumber,3)],'k','lineWidth',5);
    end
    
    % Plot 2D Map
    subplot(1,2,2)
    beamDirection= [sin(centerElevationAngle)*cos(beamAzimuthAngles(i)), ...
        sin(centerElevationAngle)*sin(beamAzimuthAngles(i))];
    data = beamDirection.*timeDomain'*soundSpeed;
    rotMatrix = [cos(-centerAzimuthAngle-pi()/2) -sin(-centerAzimuthAngle-pi()/2);
                 sin(-centerAzimuthAngle-pi()/2) cos(-centerAzimuthAngle-pi()/2)];
    dataRotated = rotMatrix*data';
    xData = dataRotated(1,:); yData = dataRotated(2,:); zData = SignalBeam(:,i);
    scatter(xData(1:length(xData)/2),yData(1:length(xData)/2), ...
        30,zData(1:length(xData)/2),'filled'); 
    hold on; %axis equal;
    
    waitbar(i/nBeams,f,['Calculating ' num2str(nBeams) ' Beams with ' num2str(nRays) ' rays']);
end

title(['Calculating ' num2str(nBeams) ' Beams with ' num2str(nRays) ' rays']);

















