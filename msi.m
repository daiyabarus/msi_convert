clear all;
close all;
clc

[Horizontal, Vertical, Optional] = msiread("80010672_1855_x_co_m45_00t.msi");
fc = Optional.frequency;


%% Normalize the azimuth to 0 degrees
Horizontal.Azimuth = mod(Horizontal.Azimuth, 360);  % Ensures values between 0 and 360 degrees
Horizontal.Azimuth = Horizontal.Azimuth - Horizontal.Azimuth(1);  % Normalize to 0 degrees

%% Pattern reconstruction and transparency adjustment
vertSlice = Vertical.Magnitude;
theta = 90 - Vertical.Elevation;  % Convert to standard spherical elevation
horizSlice = Horizontal.Magnitude;
phi = Horizontal.Azimuth;

%% Plot the radiation pattern using patternFromSlices
pOptions = PatternPlotOptions(Transparency = 0.8, AntennaOffset = [0 0 0]);
patternFromSlices(vertSlice, theta, horizSlice, phi, PatternOptions = pOptions);

%% Get the current figure handle
fig = gcf; % Get the current figure

%% Get the surface data from the current figure
hSurface = findobj(fig, 'Type', 'Surface');

if isempty(hSurface)
    error('No surface plot found in the current figure.');
end

% Get X, Y, Z data (corresponds to the radiation pattern mesh)
X = get(hSurface, 'XData');
Y = get(hSurface, 'YData');
Z = get(hSurface, 'ZData');

%% Export the 3D mesh to STL format
% You will need the surf2stl function, which you can get from MATLAB's File Exchange
surf2stl('80010672_1855_x_co_m45_00t.stl', X, Y, Z);

disp('STL file successfully generated: 80010672_1855_x_co_m45_00t.stl');
