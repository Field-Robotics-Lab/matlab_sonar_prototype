%% Beam Display

clear all
close all
clc

% distance = randn*linspace(1,60,100);
% distance = distance(randperm(length(distance)));
% alpha = randn*linspace(1,0.5236,100);
% alpha = alpha(randperm(length(alpha)));

a1 = [0 0 0];
a2 = [0.2618 0.2618 0.2618];
a3 = [0.3618 0.3618 0.3618];
a4 = [0.5 0.5 0.5];
a5 = [1 1 1];

d1 = [1 1.5 .5];
d2 = [.75 1 1.2];
d3 = [2 2.5 3];
d4 = [3 4 5];
d5 = [5 6 7];
d6 = [7 8 9];
d7 = [9 10 11];
d8 = [11 12 13];
d9 = [13 14 15];
d10 = [15 16 17];
d11 = [20 25 30];
d12 = [25 28 30];
d13 = [30 35 40];
d14 = [35 40 42];
d15 = [40 45 50];
d16 = [50 52 54];
d17 = [52 55 58];
d18 = [60 62 65];
d19 = [60 65 70];
d20 = [70 75 70];
d21 = [75 80 85];
d22 = [80 85 90];
d23 = [85 90 95];
d24 = [90 95 100];
d25 = [95 100 105];

beamMatrix = [beam(d1,a1) beam(d2,a1) beam(d3,a1) beam(d4,a1) beam(d5,a1);...
              beam(d6,a2) beam(d7,a2) beam(d8,a2) beam(d9,a2) beam(d10,a2); ...
              beam(d11,a3) beam(d12,a3) beam(d13,a3) beam(d14,a3) beam(d15,a3);...
              beam(d16,a4) beam(d17,a4) beam(d18,a4) beam(d19,a4) beam(d20,a4);...
              beam(d21,a5) beam(d22,a5) beam(d23,a5) beam(d24,a5) beam(d25,a5)];




imagesc(beamMatrix)
c = colorbar;
c.Label.String='Sound Pressure Level, [dB]';
xlabel('Horizontal Field of View')
ylabel('Vertical Field of View')