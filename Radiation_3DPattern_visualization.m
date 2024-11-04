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


%%%%%  Reshape dB Gain in 2.4GHz into a matrix %%%%%%%%%%%%%%%%%%%%%%%%
totalGainMatrixdB_2_4GHz = reshape(Gain_dB_2_4GHz, numTheta, numPhi);
%% differet normalization method
%totalGainMatrixdB_2_4GHz = totalGainMatrixdB_2_4GHz - min(min(totalGainMatrixdB_2_4GHz));

totalGainMatrix_linear_2_4GHz = 10.^(totalGainMatrixdB_2_4GHz/10);

totalGainMatrixdB_2_4GHz = totalGainMatrixdB_2_4GHz + 19;

[Phi, Theta] = meshgrid(Phi, Theta);

x = sin(Theta).*cos(Phi);
y = sin(Theta).*sin(Phi);
z = cos(Theta);

xp =  totalGainMatrixdB_2_4GHz .* x;
yp =  totalGainMatrixdB_2_4GHz .* y;
zp =  totalGainMatrixdB_2_4GHz .* z;

figure(1);

set(gcf, 'Renderer', 'painters') % Using Painters for vector graphics
%Plot3AxisAtOrigin(xp,yp,zp, 'r')
surf(xp, yp, zp, totalGainMatrix_linear_2_4GHz.*8-10, 'EdgeColor', 'none');

view(45, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
% Define the center and radius of the circles
% Define the center and radius of the circles
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 23;         % Radius of the circles

% Set axis limits to make axes touch the circle
xlim([-radius, radius]);
ylim([-radius, radius]);
zlim([-radius, radius]);

% **Arrow properties**
arrow_length = radius;
arrow_head_length = 2; % Length of the arrowhead
arrow_head_width = 1;  % Width of the arrowhead

% **Draw axes without arrowheads using plot3**
% X-axis (red)
plot3([0, arrow_length - arrow_head_length], [0, 0], [0, 0], 'g', 'LineWidth', 1);
hold on;

% Y-axis (green)
plot3([0, 0], [0, arrow_length - arrow_head_length], [0, 0], 'b', 'LineWidth', 1);

% Z-axis (blue)
plot3([0, 0], [0, 0], [0, arrow_length - arrow_head_length], 'r', 'LineWidth', 1);

% **Create solid arrowheads using 'patch'**
% X-axis arrowhead (red)
x_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
y_head = [-arrow_head_width/2, arrow_head_width/2, 0];
z_head = [0, 0, 0];
patch(x_head, y_head, z_head, 'g', 'EdgeColor', 'g');

% Y-axis arrowhead (green)
x_head = [0, 0, 0];
y_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
z_head = [-arrow_head_width/2, arrow_head_width/2, 0];
patch(x_head, y_head, z_head, 'b', 'EdgeColor', 'b');

% Z-axis arrowhead (blue)
x_head = [-arrow_head_width/2, arrow_head_width/2, 0];
y_head = [0, 0, 0];
z_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
patch(x_head, y_head, z_head, 'r', 'EdgeColor', 'r');

% Remove ticks and axis numbers
set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');
Z=get(gca,'Ztick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');
ZL=get(gca,'ZtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./30;
Yoff=diff(get(gca,'YLim'))./30;
Zoff=diff(get(gca,'ZLim'))./30;

% DRAW TICKS
%%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
for i=1:length(X)
   plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
end;
for i=1:length(Y)
   plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
end;
for i=1:length(Z)
   plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
end;

% DRAW LABELS
text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

%title('Normalized Field Pattern (2.4GHz)');% Add the colorbar below the plot
hColorbar = colorbar('southoutside'); % Place the colorbar below the plot
% Shrink the size of the colorbar (reduce its height)
cpos = hColorbar.Position; % Get the current position of the colorbar
hColorbar.Position = [cpos(1) + 0.15, cpos(2), cpos(3) - 0.37, cpos(4)];
hColorbar.FontSize = 12; % Set the colorbar font size to match the rest of the plot
hColorbar.FontName = 'Arial'; 


grid on;

% Remove the coordinate system
axis on;


% Adjust plot settings (optional)
axis equal; % Equalize axes scaling
colormap('jet'); % Choose colormap for color representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the position and size of the letters

% Define the position and size of the letter
x_pos1 = 0.05;  % X position of the letter
y_pos1 = 25;  % Y position of the letter
z_pos1 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter
font_name = 'Arial';


% Draw the letter "R" using the 'text' function
text(x_pos1, y_pos1, z_pos1,  'Y', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos3 = 0.05;  % X position of the letter
y_pos3 = 0;  % Y position of the letter
z_pos3 = 25;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos3, y_pos3, z_pos3, 'Z', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos2 = 25;  % X position of the letter
y_pos2 = 0;  % Y position of the letter
z_pos2 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos2, y_pos2, z_pos2, 'X', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);


% Define the center and radius of the circle
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 20;         % Radius of the circle

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

plot3(x, y, z, 'r', 'LineWidth', 1);  % Plot the circle
plot3(z, x, y, 'g', 'LineWidth', 1);  % Plot the circle
plot3(y, z, x, 'b', 'LineWidth', 1);  % Plot the circle

% Set plot title and labels

% Set plot labels with IEEE standard font size
xlabel('X', 'FontSize', font_size, 'FontName', font_name);
ylabel('Y', 'FontSize', font_size, 'FontName', font_name);
zlabel('Z', 'FontSize', font_size, 'FontName', font_name);

print(gcf, '-dpdf', 'Normalized_Field_Pattern_2_4GHz.pdf');

%%%%%%%%%%%%%%%%%%%%% 2.45GHz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi = unique(phiAngles);
Theta = unique(thetaAngles);
numPhi = numel(unique(phiAngles));
numTheta = numel(unique(thetaAngles));


totalGainMatrixdB_2_45GHz = reshape(Gain_dB_2_45GHz, numTheta, numPhi);
%% differet normalization method
%totalGainMatrixdB_2_4GHz = totalGainMatrixdB_2_4GHz - min(min(totalGainMatrixdB_2_4GHz));

totalGainMatrix_linear_2_45GHz = 10.^(totalGainMatrixdB_2_45GHz/10);

totalGainMatrixdB_2_45GHz = totalGainMatrixdB_2_45GHz + 19;

[Phi, Theta] = meshgrid(Phi, Theta);

x = sin(Theta).*cos(Phi);
y = sin(Theta).*sin(Phi);
z = cos(Theta);

xp =  totalGainMatrixdB_2_45GHz .* x;
yp =  totalGainMatrixdB_2_45GHz .* y;
zp =  totalGainMatrixdB_2_45GHz .* z;

figure(2);

set(gcf, 'Renderer', 'painters') % Using Painters for vector graphics
%Plot3AxisAtOrigin(xp,yp,zp, 'r')
surf(xp, yp, zp, totalGainMatrix_linear_2_45GHz.*10-11, 'EdgeColor', 'none');

view(45, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
% Define the center and radius of the circles
% Define the center and radius of the circles
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 23;         % Radius of the circles

% Set axis limits to make axes touch the circle
xlim([-radius, radius]);
ylim([-radius, radius]);
zlim([-radius, radius]);

% **Arrow properties**
arrow_length = radius;
arrow_head_length = 2; % Length of the arrowhead
arrow_head_width = 1;  % Width of the arrowhead

% **Draw axes without arrowheads using plot3**
% X-axis (red)
plot3([0, arrow_length - arrow_head_length], [0, 0], [0, 0], 'g', 'LineWidth', 1);
hold on;

% Y-axis (green)
plot3([0, 0], [0, arrow_length - arrow_head_length], [0, 0], 'b', 'LineWidth', 1);

% Z-axis (blue)
plot3([0, 0], [0, 0], [0, arrow_length - arrow_head_length], 'r', 'LineWidth', 1);

% **Create solid arrowheads using 'patch'**
% X-axis arrowhead (red)
x_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
y_head = [-arrow_head_width/2, arrow_head_width/2, 0];
z_head = [0, 0, 0];
patch(x_head, y_head, z_head, 'g', 'EdgeColor', 'g');

% Y-axis arrowhead (green)
x_head = [0, 0, 0];
y_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
z_head = [-arrow_head_width/2, arrow_head_width/2, 0];
patch(x_head, y_head, z_head, 'b', 'EdgeColor', 'b');

% Z-axis arrowhead (blue)
x_head = [-arrow_head_width/2, arrow_head_width/2, 0];
y_head = [0, 0, 0];
z_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
patch(x_head, y_head, z_head, 'r', 'EdgeColor', 'r');

% Remove ticks and axis numbers
set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');
Z=get(gca,'Ztick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');
ZL=get(gca,'ZtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./30;
Yoff=diff(get(gca,'YLim'))./30;
Zoff=diff(get(gca,'ZLim'))./30;

% DRAW TICKS
%%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
for i=1:length(X)
   plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
end;
for i=1:length(Y)
   plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
end;
for i=1:length(Z)
   plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
end;

% DRAW LABELS
text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

%title('Normalized Field Pattern (2.4GHz)');% Add the colorbar below the plot
hColorbar = colorbar('southoutside'); % Place the colorbar below the plot
% Shrink the size of the colorbar (reduce its height)
cpos = hColorbar.Position; % Get the current position of the colorbar
hColorbar.Position = [cpos(1) + 0.15, cpos(2), cpos(3) - 0.37, cpos(4)];
hColorbar.FontSize = 12; % Set the colorbar font size to match the rest of the plot
hColorbar.FontName = 'Arial'; 

grid on;

% Remove the coordinate system
axis on;


% Adjust plot settings (optional)
axis equal; % Equalize axes scaling
colormap('jet'); % Choose colormap for color representation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the position and size of the letters

% Define the position and size of the letter
x_pos1 = 0.05;  % X position of the letter
y_pos1 = 25;  % Y position of the letter
z_pos1 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter
font_name = 'Arial';


% Draw the letter "R" using the 'text' function
text(x_pos1, y_pos1, z_pos1,  'Y', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos3 = 0.05;  % X position of the letter
y_pos3 = 0;  % Y position of the letter
z_pos3 = 25;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos3, y_pos3, z_pos3, 'Z', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos2 = 25;  % X position of the letter
y_pos2 = 0;  % Y position of the letter
z_pos2 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos2, y_pos2, z_pos2, 'X', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);


% Define the center and radius of the circle
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 20;         % Radius of the circle

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

plot3(x, y, z, 'r', 'LineWidth', 1);  % Plot the circle
plot3(z, x, y, 'g', 'LineWidth', 1);  % Plot the circle
plot3(y, z, x, 'b', 'LineWidth', 1);  % Plot the circle

% Set plot title and labels

% Set plot labels with IEEE standard font size
xlabel('X', 'FontSize', font_size, 'FontName', font_name);
ylabel('Y', 'FontSize', font_size, 'FontName', font_name);
zlabel('Z', 'FontSize', font_size, 'FontName', font_name);

print(gcf, '-dpdf', 'Normalized_Field_Pattern_2_45GHz.pdf');



%%%%%%%%%%%%%%%%%%%%% 2.5GHz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi = unique(phiAngles);
Theta = unique(thetaAngles);
numPhi = numel(unique(phiAngles));
numTheta = numel(unique(thetaAngles));


totalGainMatrixdB_2_5GHz = reshape(Gain_dB_2_5GHz, numTheta, numPhi);
%% differet normalization method
%totalGainMatrixdB_2_4GHz = totalGainMatrixdB_2_4GHz - min(min(totalGainMatrixdB_2_4GHz));

totalGainMatrix_linear_2_5GHz = 10.^(totalGainMatrixdB_2_5GHz/10);


totalGainMatrixdB_2_5GHz = totalGainMatrixdB_2_5GHz + 19;

[Phi, Theta] = meshgrid(Phi, Theta);

x = sin(Theta).*cos(Phi);
y = sin(Theta).*sin(Phi);
z = cos(Theta);

xp =  totalGainMatrixdB_2_5GHz .* x;
yp =  totalGainMatrixdB_2_5GHz .* y;
zp =  totalGainMatrixdB_2_5GHz .* z;

figure(3);

set(gcf, 'Renderer', 'painters') % Using Painters for vector graphics
%Plot3AxisAtOrigin(xp,yp,zp, 'r')
surf(xp, yp, zp, totalGainMatrix_linear_2_5GHz.*10-11, 'EdgeColor', 'none');

view(45, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
% Define the center and radius of the circles
% Define the center and radius of the circles
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 23;         % Radius of the circles

% Set axis limits to make axes touch the circle
xlim([-radius, radius]);
ylim([-radius, radius]);
zlim([-radius, radius]);

% **Arrow properties**
arrow_length = radius;
arrow_head_length = 2; % Length of the arrowhead
arrow_head_width = 1;  % Width of the arrowhead

% **Draw axes without arrowheads using plot3**
% X-axis (red)
plot3([0, arrow_length - arrow_head_length], [0, 0], [0, 0], 'g', 'LineWidth', 1);
hold on;

% Y-axis (green)
plot3([0, 0], [0, arrow_length - arrow_head_length], [0, 0], 'b', 'LineWidth', 1);

% Z-axis (blue)
plot3([0, 0], [0, 0], [0, arrow_length - arrow_head_length], 'r', 'LineWidth', 1);

% **Create solid arrowheads using 'patch'**
% X-axis arrowhead (red)
x_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
y_head = [-arrow_head_width/2, arrow_head_width/2, 0];
z_head = [0, 0, 0];
patch(x_head, y_head, z_head, 'g', 'EdgeColor', 'g');

% Y-axis arrowhead (green)
x_head = [0, 0, 0];
y_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
z_head = [-arrow_head_width/2, arrow_head_width/2, 0];
patch(x_head, y_head, z_head, 'b', 'EdgeColor', 'b');

% Z-axis arrowhead (blue)
x_head = [-arrow_head_width/2, arrow_head_width/2, 0];
y_head = [0, 0, 0];
z_head = [arrow_length - arrow_head_length, arrow_length - arrow_head_length, arrow_length];
patch(x_head, y_head, z_head, 'r', 'EdgeColor', 'r');

% Remove ticks and axis numbers
set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');
Z=get(gca,'Ztick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');
ZL=get(gca,'ZtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./30;
Yoff=diff(get(gca,'YLim'))./30;
Zoff=diff(get(gca,'ZLim'))./30;

% DRAW TICKS
%%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
for i=1:length(X)
   plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
end;
for i=1:length(Y)
   plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
end;
for i=1:length(Z)
   plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
end;

% DRAW LABELS
text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

%title('Normalized Field Pattern (2.4GHz)');% Add the colorbar below the plot
hColorbar = colorbar('southoutside'); % Place the colorbar below the plot
% Shrink the size of the colorbar (reduce its height)
cpos = hColorbar.Position; % Get the current position of the colorbar
hColorbar.Position = [cpos(1) + 0.15, cpos(2), cpos(3) - 0.37, cpos(4)];
hColorbar.FontSize = 12; % Set the colorbar font size to match the rest of the plot
hColorbar.FontName = 'Arial'; 

grid on;

% Remove the coordinate system
axis on;


% Adjust plot settings (optional)
axis equal; % Equalize axes scaling
colormap('jet'); % Choose colormap for color representation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the position and size of the letters

% Define the position and size of the letter
x_pos1 = 0.05;  % X position of the letter
y_pos1 = 25;  % Y position of the letter
z_pos1 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter
font_name = 'Arial';


% Draw the letter "R" using the 'text' function
text(x_pos1, y_pos1, z_pos1,  'Y', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos3 = 0.05;  % X position of the letter
y_pos3 = 0;  % Y position of the letter
z_pos3 = 25;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos3, y_pos3, z_pos3, 'Z', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);

% Define the position and size of the letter
x_pos2 = 25;  % X position of the letter
y_pos2 = 0;  % Y position of the letter
z_pos2 = 0.05;  % Z position of the letter
font_size = 15;  % Font size of the letter

% Draw the letter "R" using the 'text' function
text(x_pos2, y_pos2, z_pos2, 'X', 'Interpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);


% Define the center and radius of the circle
center = [0, 0, 0];  % Center coordinates [x, y, z]
radius = 20;         % Radius of the circle

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

plot3(x, y, z, 'r', 'LineWidth', 1);  % Plot the circle
plot3(z, x, y, 'g', 'LineWidth', 1);  % Plot the circle
plot3(y, z, x, 'b', 'LineWidth', 1);  % Plot the circle

% Set plot title and labels

% Set plot labels with IEEE standard font size
xlabel('X', 'FontSize', font_size, 'FontName', font_name);
ylabel('Y', 'FontSize', font_size, 'FontName', font_name);
zlabel('Z', 'FontSize', font_size, 'FontName', font_name);

print(gcf, '-dpdf', 'Normalized_Field_Pattern_2_5GHz.pdf');
