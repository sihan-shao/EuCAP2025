clc;clear;
Phi = importdata("2023 AUT Cantenna_Realized Gain dBi - Phi component - 3GHz.txt");
Phi = Phi.data;
Theta = importdata("2023 AUT Cantenna_Realized Gain dBi - Theta component - 3GHz.txt");
Theta = Theta.data;
% Extract the angles and gains
phiAngles = Phi(:, 1);
thetaAngles = Phi(:, 2);
thetaGains = Phi(:, 3);
phiGains = Theta(:, 3);
% Convert the gains from dB to linear scale
thetaGainsLinear = 10.^(thetaGains / 10);
phiGainsLinear = 10.^(phiGains / 10);
% Calculate the total gain in linear scale
totalGainLinear = thetaGainsLinear + phiGainsLinear;
% Convert the total gain back to dB scale
totalGain = 10 * log10(totalGainLinear);


% Determine the number of unique phi and theta angles
Phi = unique(phiAngles);
Theta = unique(thetaAngles);
numPhi = numel(unique(phiAngles));
numTheta = numel(unique(thetaAngles));

% Reshape totalGainLinear into a matrix
totalGainMatrix = reshape(totalGainLinear, numTheta, numPhi);
totalGainMatrixdB = reshape(totalGain, numTheta, numPhi);
totalGainMatrixdB = totalGainMatrixdB - min(min(totalGainMatrixdB));

% % Create a meshgrid for 3D plotting
[Phi, Theta] = meshgrid(Phi, Theta);
[x,y,z]=sph2cart(Phi,pi/2-Theta,totalGainMatrix);
figure(1);
surf(x,y,z);
shading interp;
title('3D Radiation Pattern for Total Gain (Absolute Value)');
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap("jet");
colorbar;

[x,y,z]=sph2cart(Phi,pi/2-Theta,totalGainMatrixdB);
figure(2);
surf(x,y,z);
shading interp;
title('3D Radiation Pattern for Total Gain (dB)');
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap("jet");
colorbar;


% Reshape totalGain into a matrix
totalGainMatrix = reshape(totalGain, numTheta, numPhi);
Phi = unique(phiAngles);
Theta = unique(thetaAngles);

% Find the indices where Theta equals -90/90
theta90Indices1 = find(abs(Theta + pi/2) < 0.0001);
theta90Indices2 = find(abs(Theta - pi/2) < 0.0001);
GaindB1 = [totalGainMatrix(theta90Indices1,:) totalGainMatrix(theta90Indices2,:)];
rho1 = linspace(-pi, pi, length(GaindB1));
figure(3);
polarplot(rho1,GaindB1);
title('Antenna Radiation Pattern (H-plane) for \theta=90^o (dB)')
rmin = min(GaindB1);
rmax = max(GaindB1);
rlim([rmin rmax])


phi0Indices = find(abs(Phi) < 0.0001);
GaindB2 = totalGainMatrix(:,phi0Indices);
rho = Theta;
figure(4)
polarplot(rho, GaindB2);
hold on
title('Antenna Radiation Pattern (XZ-plane) for \phi=0^o (dB)')
rmin = min(GaindB2);
rmax = max(GaindB2);
rlim([rmin rmax])

% Calculate the -3dB point
dB3 = max(GaindB2) - 3;
% Find the indices where the gain is close to the -3dB point
dB3Indices = find(abs(GaindB2 - dB3) < 0.1);
% Extract the corresponding Theta angles
thetaDB3 = rho(dB3Indices);
% Draw lines at the -3dB points
polarplot([thetaDB3(1) thetaDB3(1)], rlim, 'r--');
polarplot([thetaDB3(2) thetaDB3(2)], rlim, 'r--');

% Calculate the bandwidth
HPBW = abs(thetaDB3(1) - thetaDB3(2));

% Draw the bandwidth value on the figure
HPBWDegrees = rad2deg(HPBW);
text(0, 0, ['HPBW: ' num2str(HPBWDegrees) '^o'], 'HorizontalAlignment', 'center');


figure(5);
plot(rho, GaindB2);
hold on
title('Antenna Radiation Pattern (XZ-plane) for \phi=0^o (dB)')
ymin = min(GaindB2);
ymax = max(GaindB2);
ylim([ymin ymax]);
xlim([-pi pi]);
plot([thetaDB3(1) thetaDB3(1)], ylim, 'r--')
plot([thetaDB3(2) thetaDB3(2)], ylim, 'r--')
dB3 = ymax - 3;
plot(xlim, [dB3 dB3], 'b--');
xlabel('Theta (rad.)');
ylabel('Realized Gain (dB)');
text(0, 0, ['HPBW: ' num2str(HPBWDegrees) '^o'], 'HorizontalAlignment', 'center');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0Indices = find(abs(Phi-pi/2) < 0.0001);
GaindB3 = totalGainMatrix(:,phi0Indices);
rho = Theta;
figure(6)
polarplot(rho, GaindB3);
hold on
title('Antenna Radiation Pattern (YZ-plane) for \phi=90^o (dB)')
rmin = min(GaindB3);
rmax = max(GaindB3);
rlim([rmin rmax])

% Calculate the -3dB point
dB3 = max(GaindB3) - 3;
% Find the indices where the gain is close to the -3dB point
dB3Indices = find(abs(GaindB3 - dB3) < 0.2);
% Extract the corresponding Theta angles
thetaDB3 = rho(dB3Indices);
% Draw lines at the -3dB points
polarplot([thetaDB3(1) thetaDB3(1)], rlim, 'r--');
polarplot([thetaDB3(2) thetaDB3(2)], rlim, 'r--');

% Calculate the bandwidth
HPBW = abs(thetaDB3(1) - thetaDB3(2));

% Draw the bandwidth value on the figure
HPBWDegrees = rad2deg(HPBW);
text(0, 0, ['HPBW: ' num2str(HPBWDegrees) '^o'], 'HorizontalAlignment', 'center');


figure(7);
plot(rho, GaindB3);
hold on
title('Antenna Radiation Pattern (YZ-plane) for \phi=90^o (dB)')
ymin = min(GaindB3);
ymax = max(GaindB3);
ylim([ymin ymax]);
xlim([-pi pi]);
plot([thetaDB3(1) thetaDB3(1)], ylim, 'r--')
plot([thetaDB3(2) thetaDB3(2)], ylim, 'r--')
dB3 = ymax - 3;
plot(xlim, [dB3 dB3], 'b--');
xlabel('Theta (rad.)');
ylabel('Realized Gain (dB)');
text(0, 0, ['HPBW: ' num2str(HPBWDegrees) '^o'], 'HorizontalAlignment', 'center');

