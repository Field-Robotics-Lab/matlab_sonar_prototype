%% Single Beam Sonar
%% Sonar Point-Scatter Model
% Contributors: Andreina Rascon, Derek Olson, Woeng-Sug Choi

clear;
close all;
clc;

%% Sonar Parameters
% Input Parameters
% The BlueView P900-45 is used as an example sonar for the purposes of
% demonstrating the point-scatter model

sonarFreq = 900E3; % Sonar frquency
bandwidth = 29.5e4; % [Hz]
freqResolution = 100e2;

nBeams = 1;
beam.azimuthAngle = 0.0; % Beam at center line in azimuth direction
beam.elevationAngle = 0.0175; % Beam looking down in elevation direction
beam.azimuthAngleWidth = 0.1; % radians
beam.elevationAngleWidth = 0.1; % radians

ray.nAzimuthRays = 4;
ray.nElevationRays = 3;

%% Input Parameters
soundSpeed = 1500.0; % [m/s]

% Target information (which should be obtained from ray gazebo plugin)
% 3 Targets, no incident angle

ray.distance = [15,5,10;2,10,10;15,15,15;4,2,3];
ray.alpha = [0,0,0;0,0,0;0,0,0;0,0,0]; % no incident angle
ray.azimuthAngles = beam.azimuthAngle + linspace(-beam.azimuthAngleWidth/2,beam.azimuthAngleWidth/2,ray.nAzimuthRays);
ray.elevationAngles = beam.elevationAngle + linspace(-beam.elevationAngleWidth/2,beam.elevationAngleWidth/2,ray.nElevationRays);
ray.azimuthAngleWidth = beam.azimuthAngleWidth/(ray.nAzimuthRays-1);
ray.elevationAngleWidth = beam.elevationAngleWidth/(ray.nElevationRays-1);

% Surface Reflectivity
mu = 10e-4;

% Frequency spectrum
fmin = sonarFreq - bandwidth/2*4;
fmax = sonarFreq + bandwidth/2*4;

% Sampling periods
max_T = max(max(ray.distance)*2/soundSpeed);
delta_f = 1/max_T;
delta_t = 1/(fmax - fmin);

nFreq = round((fmax - fmin) / delta_f);
freq = linspace(fmin,fmax,nFreq);
time = linspace(0,max_T,nFreq);

% Transmit spectrum, frequency domain
S_f = 1e11*exp(-(freq-sonarFreq).^2.*pi^2./bandwidth^2);

% define wave vector
kw = 2*pi*freq/soundSpeed;
        
% add attenuation
absorption = 0.0354; % [dB/m]
attenuation = absorption*log(10)/20;
K = kw + 1j*attenuation;

%% Point Scattering model
% Calculate echo level using the point scatter model

P_ray_f = zeros(ray.nAzimuthRays,ray.nElevationRays,nFreq); % pre-allocate P(f) for a beam
P_ray_t = zeros(ray.nAzimuthRays,ray.nElevationRays,nFreq); % pre-allocate P(t) for a beam

% --- Beam loop --- %
% Azimuth ray Loop
for k = 1:ray.nAzimuthRays
    % Elevation ray Loop
    for i = 1:ray.nElevationRays
        xi_z = rand;   % generate a random number, (Gaussian noise)
        xi_y = rand;   % generate another random number, (Gaussian noise)
        
        alpha = ray.alpha(k,i);   % angle between ray vector and object normal vector, [rad]
        distance = ray.distance(k,i);
        azimuthBeamPattern = (sinc(pi()*0.884/ray.azimuthAngleWidth*sin(ray.azimuthAngles(k)))).^2;
        elevationBeamPattern = (sinc(pi()*0.884/ray.elevationAngleWidth*sin(ray.elevationAngles(i)))).^2;
        amplitude = ((xi_z+1i*xi_y)/sqrt(2))*(sqrt(mu*cos(alpha)^2*distance^2*ray.azimuthAngleWidth*ray.elevationAngleWidth))*azimuthBeamPattern*elevationBeamPattern;
        
        % Summation of Echo returned from a signal (frequency domain)
        for m = 1:nFreq
            P_ray_f(k,i,m) = P_ray_f(k,i,m) + S_f(m)*amplitude*exp(-1i*K(m)*distance*2)/(distance^2);
        end
    end % end of elevation ray loop
end % end of azimuth ray loop

% Summation for a beam
P_beam = squeeze(sum(sum(P_ray_f,2),1));
P_beam_t = ifft(P_beam);

%% Plotting

% Plot inverse fast fourier transform 
% figure (1)
% figure;
% plot(time,P_beam_t)
% xlabel('Time, [s]')
% ylabel('Pressure, [Pa]')


% Plot Sound Pressure Level of Echo Level
SPL = 20*log(abs(P_beam_t)); % sound pressure level, [dB]

figure;
plot(time,SPL)
xlabel('Time, [s]')
ylabel('Sound Pressure Level, [dB]')
grid on; hold on;
% axis([0 3e-3 -400 400])
set(0, 'defaulttextfontsize',20)
set(0, 'defaultaxesfontsize',10)
set(gcf,'color','w')
% export_fig SingleBeam_SPL.png -png -r300 -painters



function y = sinc(x)
     y = abs(sin(x)./x);
     y(x==0) = 1;
end

