clear, clc, close all

addpath('function')

% Set parameters
BS.loc  = [-25; 25; 15];
BS.ori  = [315; 12; 0];
BS.H    = 8;
BS.V    = 8;

UE.ori  = [0; 0; 0];
UE.H    = 1;
UE.V    = 1;
UE.h    = 1.5;

% Load UE locations
load('location_xs.mat')                 % location (x_UE, y_UE) [m], dimensions:(# of UE) x 2
location(4e5 + 1 : end, :) = [];
location(:, 3)  = UE.h;                 % location (x_UE, y_UE, h_UE) [m], dimensions: (# of UE) x 3
len_channel     = size(location, 1);
location        = transpose(location);  % Dimensions: 3 x (# of UE)

% Generate channels using Saleh-Valenzuela (SV) channel model
ch.L           = 20;               % The number of multi-paths
ch.sigma_alpha = 1;                % Std. dev. of complex gain for channels
ch.sigma_L     = 0;                % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
ch.radiation_pattern = false;      % Radiation power pattern: ture - The radiation power pattern based on 3GPP TR 38.901, false - Isotropic power pattern

channel = zeros(len_channel, BS.H*BS.V);
for i = 1 : size(location, 2)
    UE.loc = location(:, i);
    channel(i, :) = SV_channel(BS, UE, ch);
end

% Save the channel
if ~isfolder('channel')
    mkdir('channel')
    addpath('channel')
end
save(['channel/H', num2str(BS.H), 'V', num2str(BS.V), '_sigmaL_', num2str(ch.sigma_L)], 'channel')