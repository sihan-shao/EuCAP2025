clc
clear

% Specify the file path
H_file_path = '2023 AUT Cantenna_Realized Gain dBi - Phi component - 3GHz.txt';

% Import the data
H_data = importdata(H_file_path);

% Access each column as a 1D array
H_column1 = H_data.data(:, 1);
H_column2 = H_data.data(:, 2);
H_gain_dB = H_data.data(:, 3);

% Display the imported data
% disp(H_column1);
% disp(H_column2);
% disp(H_column3);

% Specify the file path
V_file_path = '2023 AUT Cantenna_Realized Gain dBi - Theta component - 3GHz.txt';

% Import the data
V_data = importdata(V_file_path);

% Access each column as a 1D array
V_column1 = V_data.data(:, 1);
V_column2 = V_data.data(:, 2);
V_gain_dB = V_data.data(:, 3);

% Display the imported data
% disp(V_column1);
% disp(V_column2);
% disp(V_column3);

% Assuming you have the gain values stored in a variable named 'gain_dB'
% Calculate the corresponding absolute gain in linear scale
H_gain_linear = 10.^(H_gain_dB/10);
V_gain_linear = 10.^(V_gain_dB/10);
T_gain_linear = H_gain_linear + V_gain_linear;
% Display the absolute gain values
% disp(H_gain_linear);
% disp(V_gain_linear);
% disp(T_gain_linear);

T_gain_dB = 10*log10(T_gain_linear);


% Example 1D gain array with 7260 elements
%T_gain_dB = [gain1, gain2, gain3, ..., gain7260];

% Number of phi values
num_phi = 60; % Assuming you have 60 unique phi values

% Reshape the 1D gain array into a 2D array
T_gain_dB_2d = reshape(T_gain_dB, 121, []);
% disp(T_gain_dB_2d);
% Reshape the 2D array to have each phi value in a separate column
% gain_2d = reshape(gain_2d, [], num_phi);
T_gain_linear_2d = 10.^(T_gain_dB_2d/10);

theta_rad = H_column2(1:121);
phi_rad = H_column1(1:121:end);
gain = T_gain_linear_2d ./ max(T_gain_linear_2d) ;
size(gain)
disp(max(T_gain_linear_2d(:)));
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-----------------------------Plotting------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% % Measurement data (replace with your actual data)
% theta = H_column1; % Theta values
% phi = H_column2; % Phi values
% gain = T_gain_dB_2d; % Gain values
% 
% % Convert theta and phi to radians
% theta_rad = deg2rad(theta);
% phi_rad = deg2rad(phi);
% 
% % Convert gain to dB
% gain_dB = 10*log10(gain);
% 


% Convert theta, phi, and gain to Cartesian coordinates
[theta_grid, phi_grid] = meshgrid(theta_rad, phi_rad);
x = sin(theta_grid).*cos(phi_grid) ; % Size: 121x60
y = sin(theta_grid).*sin(phi_grid) ; % Size: 121x60
z = cos(theta_grid) ; % Size: 121x60

xp =  sqrt(gain) .* x';
yp =  sqrt(gain) .* y';
zp =  sqrt(gain) .* z';

%interpolation






% Plot 3D gain pattern
figure;
surf(xp, yp, zp, gain, 'EdgeColor', 'none');
title('Normalized Field Pattern');
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar; % Add colorbar to show gain values

% Hide the grid lines
grid on;

% Remove the coordinate system
axis on;



% Add vector arrows on the tips of the axes

% arrrow([-1.5 0 0],[1.5 0 0],'color','green','stemWidth',0.01,'facealpha',0.5)
% arrrow([0 -1.5 0],[0 1.5 0],'color','blue','stemWidth',0.01,'facealpha',0.5)
% arrrow([0 0 -1.5],[0 0 1.5],'color','red','stemWidth',0.01,'facealpha',0.5)


% Adjust plot settings (optional)
axis equal; % Equalize axes scaling
colormap('jet'); % Choose colormap for color representation


% Define the position and size of the letter
x_pos1 = 0.05;  % X position of the letter
y_pos1 = 1.5;  % Y position of the letter
z_pos1 = 0.05;  % Z position of the letter
font_size = 20;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos1, y_pos1, z_pos1, '$\mathcal{Y}$', 'Interpreter', 'latex', 'FontSize', font_size);

% Define the position and size of the letter
x_pos3 = 0.05;  % X position of the letter
y_pos3 = 0;  % Y position of the letter
z_pos3 = 1.5;  % Z position of the letter
font_size = 20;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos3, y_pos3, z_pos3, '$\mathcal{Z}$', 'Interpreter', 'latex', 'FontSize', font_size);

% Define the position and size of the letter
x_pos2 = 1.5;  % X position of the letter
y_pos2 = 0;  % Y position of the letter
z_pos2 = 0.05;  % Z position of the letter
font_size = 20;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos2, y_pos2, z_pos2, '$\mathcal{X}$', 'Interpreter', 'latex', 'FontSize', font_size);


% Define the center and radius of the circle
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 1;         % Radius of the circle

% Define the angle range for the circle (in radians)
theta = linspace(0, 2*pi, 100);

% Compute the x, y, and z coordinates of the circle points
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);
z = center(3) * ones(size(theta));  % Constant z-coordinate (for a circle in xy-plane)

% Plot the circle in 3D

hold on;
grid on;
axis equal;

plot3(x, y, z, 'r', 'LineWidth', 3);  % Plot the circle
plot3(z, x, y, 'g', 'LineWidth', 3);  % Plot the circle
plot3(y, z, x, 'b', 'LineWidth', 3);  % Plot the circle

% Set plot title and labels

xlabel('X');
ylabel('Y');
zlabel('Z');

% axis tight;
% grid on;
% % Set the number of frames for the animation
% num_frames = 150; % Adjust as desired
% 
% % Preallocate the frames array
% frames(num_frames) = struct('cdata', [], 'colormap', []);
% 
% % Create the rotating animation
% for frame = 1:num_frames
%     % Rotate the plot
%     view(frame * 360 / num_frames, 20); % Adjust the viewing angle as desired
%     
%     % Capture the current frame
%     frames(frame) = getframe(gcf);
% end
% 
% % Create a GIF file from the captured frames
% filename = 'rotation_animation.gif'; % Specify the file name
% for frame = 1:num_frames
%     % Convert the frame to an indexed image
%     [imind, cm] = rgb2ind(frames(frame).cdata, 256);
%     
%     % Write the indexed image to the GIF file
%     if frame == 1
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
%     else
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%     end
% end
% 
% disp('GIF animation generated successfully!');


figure;
% Define the theta values
pptheta = theta_rad; % Range of theta values (0 to 2*pi) with 100 points
pp2theta = theta_rad;
% Define the radius values
ppradius = 10*log10(T_gain_linear(1:121));; % Example radius values, you can replace it with your own data
pp2radius = 10*log10(T_gain_linear(7140:7260)) ;
% Create the polar plot
polarplot(pptheta, ppradius, 'LineWidth', 2, 'DisplayName', 'Phi = 0'); % Adjust line width as desired
rlim([-35 10])
hold on;
polarplot(pp2theta, pp2radius, '--' , 'LineWidth', 2, 'DisplayName', 'Phi = 180','Color','cyan');; % Adjust line width as desired
rlim([-30 10])
legend('Location', 'southwest');
% Set the plot title (optional)
title('Realized Gain (dB): XZ-plane: \Phi = 0^\circ , θ variable');
grid('minor');
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
figure;
plot(rad2deg(pptheta), ppradius, 'LineWidth', 2, 'DisplayName', 'Phi = 0');
grid on;
xlim([-180 180]);
xticks(-180:30:180);
grid('minor');
title('XZ-plane: \Phi = 0^\circ , θ variable');
xlabel('θ (deg)')
ylabel('Realized Gain (dB)')
line([-32 -32], [-25 10], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line([+33 +33], [-25 10], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line( [-180 180],[5.2 5.2], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line( [-180 180],[8.2 8.2],  'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
figure;
% Define the theta values
pp3theta = theta_rad; % Range of theta values (0 to 2*pi) with 100 points
% Define the radius values
pp3radius = 10*log10(T_gain_linear(3631:3751)); % Example radius values, you can replace it with your own data
% Create the polar plot
polarplot(pp3theta, pp3radius, 'LineWidth', 2, 'DisplayName', 'Phi = 0', 'Color','green'); % Adjust line width as desired
rlim([-30 10]);
title('Realized Gain (dB): YZ-plane: \Phi = 90^\circ , θ variable');
grid('minor');
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
figure;
plot(rad2deg(pp3theta), pp3radius, 'LineWidth', 2, 'DisplayName', 'Phi = 0', 'Color','green')
grid on;
title('YZ-plane: \Phi = 90^\circ , θ variable');
xlim([-180 180]);
xticks(-180:30:180);

xlabel('θ (deg)')
ylabel('Realized Gain (dB)')
grid('minor');
line([-32 -32], [-25 10], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line([+34 +34], [-25 10], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line( [-180 180],[5.2 5.2], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);
line( [-180 180],[8.2 8.2],  'Color', 'black', 'LineStyle', '--', 'LineWidth', 1);


% Adjust plot settings (optional)
ax = gca; % Get the current axes handle
ax.FontSize = 12; % Adjust font size as desired
