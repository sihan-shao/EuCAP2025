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

% % Calculate the gain in dB scale from two componenets
% Gain_Theta_linear_2_45GHz = data_2_45GHz.GainTheta_Lin;
% Gain_Phi_linear_2_45GHz = data_2_45GHz.GainPhi_Lin;
% Gain_Linear_2_45GHz = Gain_Theta_linear_2_45GHz + Gain_Phi_linear_2_45GHz;
% Gain_dB_calculated_2_45GHz = 10.*log10(Gain_Linear_2_45GHz);
% %Gain_dB_2_45GHz = 10.*log10(Gain_Linear_2_45GHz);
% 
% figure (1);
% 
% plot(abs(Gain_dB_calculated_2_45GHz - Gain_dB_2_45GHz));

% Determine the number of unique phi and theta angles
Phi = unique(phiAngles);
Theta = unique(thetaAngles);
numPhi = numel(unique(phiAngles));
numTheta = numel(unique(thetaAngles));

% Reshape totalGain into a matrix
totalGainMatrix_2_4GHz = reshape(Gain_dB_2_4GHz, numTheta, numPhi);
totalGainMatrix_2_45GHz = reshape(Gain_dB_2_45GHz, numTheta, numPhi);
totalGainMatrix_2_5GHz = reshape(Gain_dB_2_5GHz, numTheta, numPhi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0Indices = find(abs(Phi) < 0.0001);
GaindB_Phi_0_2_45GHz = totalGainMatrix_2_45GHz(:,phi0Indices);
rho = Theta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi90Indices = find(abs(Phi-pi/2) < 0.0001);
GaindB_Phi_90_2_45GHz = totalGainMatrix_2_45GHz(:,phi90Indices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (1);

subplot(1, 2, 2);
p1 = polarplot(rho,GaindB_Phi_0_2_45GHz, 'b-', 'LineWidth',2);
hold on;
rlim([-25 0]);
set(gca,'ThetaZeroLocation','top')

p1 = polarplot(rho,GaindB_Phi_90_2_45GHz, 'r-', 'LineWidth',2);
hold on;
rlim([-25 0]);
set(gca,'ThetaZeroLocation','top')

% Set font properties for subplot 1
ax2 = gca;
ax2.FontSize = 16;
ax2.FontName = 'Arial';


% Global settings
%title('Antenna Radiation Pattern for \phi=0^o and \phi=90^o (dB)')
%legend([p1, p2], {'\phi=0^o', '\phi=90^o'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and process the data
filename = 'Phi_0.txt';
data_phi_0 = readtable(filename, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, ...
    'ReadVariableNames', false, 'VariableNamingRule', 'preserve');

% Clean up the data
data_phi_0.Var1 = [];
data_phi_0.Var9 = str2double(data_phi_0.Var9);
data_phi_0.Properties.VariableNames = {'Theta', 'Phi', 'Abs_Grlz_dBi', 'Abs_Theta_dBi', 'Phase_Theta_deg', 'Abs_Phi_dBi', 'Phase_Phi_deg', 'Ax_Ratio_dB'};


% Extract the angles and gain
Gain_dB = data_phi_0.Abs_Grlz_dBi;          % Gain in dBi

% Convert angles to radians for polar plotting
theta = linspace(0, 2*pi, length(Gain_dB(:, 1)));

subplot(1, 2, 1);
p1 = polarplot(theta, Gain_dB, 'b-', 'LineWidth',2);
hold on;
rlim([-25 0]);
set(gca,'ThetaZeroLocation','top')

% Load and process the data (as in your original code)
filename = 'Phi_90.txt';
data_phi_90 = readtable(filename, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, ...
    'ReadVariableNames', false, 'VariableNamingRule', 'preserve');

% Clean up the data
data_phi_90.Var1 = [];
data_phi_90.Var9 = str2double(data_phi_90.Var9);
data_phi_90.Properties.VariableNames = {'Theta', 'Phi', 'Abs_Grlz_dBi', 'Abs_Theta_dBi', 'Phase_Theta_deg', 'Abs_Phi_dBi', 'Phase_Phi_deg', 'Ax_Ratio_dB'};


% Extract the angles and gain
Gain_dB = data_phi_90.Abs_Grlz_dBi;          % Gain in dBi

% Convert angles to radians for polar plotting
theta = linspace(0, 2*pi, length(Gain_dB(:, 1)));


p2 = polarplot(theta, Gain_dB, 'r-', 'LineWidth',2);
hold on;
rlim([-25 0]);
set(gca,'ThetaZeroLocation','top')

% Set font properties for subplot 2
ax1 = gca;
ax1.FontSize = 16;
ax1.FontName = 'Arial';

% Global settings
%title('Antenna Radiation Pattern for \phi=0^o and \phi=90^o (dB)')
legend([p1, p2], {'\phi=0^o', '\phi=90^o'}, 'Location', 'best');
% Adjust figure properties
set(gcf, 'Renderer', 'painters');

% Save the figure as a vectorized PDF
print(gcf, '-dpdf', 'combined_polar_plots.pdf', '-bestfit', '-r300', '-vector')