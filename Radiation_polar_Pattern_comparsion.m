clc;clear;

data = readtable('Ali_Sihan_Starlab_Calibrated.txt');

% Extract the different frequency points
data_2_4GHz = data(data.Frequency == 2.4*1e9, :);
data_2_45GHz = data(data.Frequency == 2.45*1e9, :);
data_2_5GHz = data(data.Frequency == 2.5*1e9, :);

% Extract the angles
phiAngles = data_2_4GHz.Phi;
thetaAngles = data_2_4GHz.Theta;

% Extract the gain in dB scale
Gain_dB_2_4GHz = data_2_45GHz.Gain_DB;% data_2_4GHz.Gain_DB;
Gain_dB_2_45GHz = data_2_45GHz.Gain_DB;
Gain_dB_2_5GHz = data_2_5GHz.Gain_DB;

% Determine the number of unique phi and theta angles
Phi = unique(phiAngles);
Theta = unique(thetaAngles);
numPhi = numel(unique(phiAngles));
numTheta = numel(unique(thetaAngles));

% Reshape totalGain into a matrix
totalGainMatrix_2_4GHz = reshape(Gain_dB_2_4GHz, numTheta, numPhi);
totalGainMatrix_2_45GHz = reshape(Gain_dB_2_45GHz, numTheta, numPhi);
totalGainMatrix_2_5GHz = reshape(Gain_dB_2_5GHz, numTheta, numPhi);



% Find the indices where Theta equals -90/90
theta90Indices1 = find(abs(Theta + pi/2) < 0.00001);
theta90Indices2 = find(abs(Theta - pi/2) < 0.00001);

GaindB_theta_90_2_4GHz = [totalGainMatrix_2_4GHz(theta90Indices1,:) totalGainMatrix_2_4GHz(theta90Indices2,:)];
rho1 = linspace(-pi, pi, length(GaindB_theta_90_2_4GHz));
GaindB_theta_90_2_45GHz = [totalGainMatrix_2_45GHz(theta90Indices1,:) totalGainMatrix_2_45GHz(theta90Indices2,:)];
GaindB_theta_90_2_5GHz = [totalGainMatrix_2_5GHz(theta90Indices1,:) totalGainMatrix_2_5GHz(theta90Indices2,:)];


figure(1);
p1 = polarplot(rho1,GaindB_theta_90_2_4GHz, 'r-', 'LineWidth',3);
hold on;
rlim([-18.8656 -1.0310]);

p2 = polarplot(rho1,GaindB_theta_90_2_45GHz, 'g--', 'LineWidth',3);
hold on;
rlim([-18.8656 -1.0310]);

p3 = polarplot(rho1,GaindB_theta_90_2_5GHz, 'b-.', 'LineWidth',3);
hold on;
rlim([-18.8656 -1.0310]);

% Global settings
title('Antenna Radiation Pattern (YZ-plane) for \theta=90^o (dB)');
legend([p1, p2, p3], {'2.4 GHz', '2.45 GHz', '2.5 GHz'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0Indices = find(abs(Phi) < 0.0001);
GaindB_Phi_0_2_4GHz = totalGainMatrix_2_4GHz(:,phi0Indices);
GaindB_Phi_0_2_45GHz = totalGainMatrix_2_45GHz(:,phi0Indices);
GaindB_Phi_0_2_5GHz = totalGainMatrix_2_5GHz(:,phi0Indices);

rho = Theta;
figure(2);
p1 = polarplot(rho,GaindB_Phi_0_2_4GHz, 'r-', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_0_2_4GHz) max(GaindB_Phi_0_2_4GHz)]);

p2 = polarplot(rho,GaindB_Phi_0_2_45GHz, 'g--', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_0_2_45GHz) max(GaindB_Phi_0_2_45GHz)]);

p3 = polarplot(rho,GaindB_Phi_0_2_5GHz, 'b-.', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_0_2_5GHz) max(GaindB_Phi_0_2_5GHz)]);

% Global settings
title('Antenna Radiation Pattern (XZ-plane) for \phi=0^o (dB)')
legend([p1, p2, p3], {'2.4 GHz', '2.45 GHz', '2.5 GHz'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi90Indices = find(abs(Phi-pi/2) < 0.0001);
GaindB_Phi_90_2_4GHz = totalGainMatrix_2_4GHz(:,phi90Indices);
GaindB_Phi_90_2_45GHz = totalGainMatrix_2_45GHz(:,phi90Indices);
GaindB_Phi_90_2_5GHz = totalGainMatrix_2_5GHz(:,phi90Indices);

rho = Theta;
figure(3);
p1 = polarplot(rho,GaindB_Phi_90_2_4GHz, 'r-', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_90_2_4GHz) max(GaindB_Phi_90_2_4GHz)]);

p2 = polarplot(rho,GaindB_Phi_90_2_45GHz, 'g--', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_90_2_45GHz) max(GaindB_Phi_90_2_45GHz)]);

p3 = polarplot(rho,GaindB_Phi_90_2_5GHz, 'b-.', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_90_2_5GHz) max(GaindB_Phi_90_2_5GHz)]);

% Global settings
title('Antenna Radiation Pattern (YZ-plane) for \phi=90^o (dB)')
legend([p1, p2, p3], {'2.4 GHz', '2.45 GHz', '2.5 GHz'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
p1 = polarplot(rho,GaindB_Phi_0_2_45GHz, 'b-', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_0_2_45GHz) max(GaindB_Phi_0_2_45GHz)]);


p2 = polarplot(rho,GaindB_Phi_90_2_45GHz, 'r-', 'LineWidth',3);
hold on;
rlim([min(GaindB_Phi_90_2_45GHz) max(GaindB_Phi_90_2_45GHz)]);

% Global settings
title('Antenna Radiation Pattern for \phi=0^o and \phi=90^o (dB)')
legend([p1, p2], {'\phi=0^o', '\phi=90^o'});


