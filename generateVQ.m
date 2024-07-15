addpath('function')

% Set parameters
BS.loc  = [-25; 25; 15];   % The location of the base station, (x_BS, y_BS, h_BS) [m]
BS.ori  = [315; 12; 0];    % The orientation(=bore-sight) of the base station, (the bearing angle, the down-tilt angle, the slant angle) [deg]
BS.H    = 8;               % The number of horizontal elements in the base station
BS.V    = 8;               % The number of vertical elements in the base station
UE.h    = 1.5;             % The height of UE [m]
Q       = 256;              % Codebook size

% Load UE locations
load('location_xs.mat')    % location (x_UE, y_UE) [m], dimensions:(# of UE) x 2
location(:, 3)  = UE.h;     % location (x_UE, y_UE, h_UE) [m], dimensions: (# of UE) x 3

% Remove duplicate location coordinates to reduce computational load in function of "Codebook_VQ".
location        = round(location, 3);

% Generate the perfect location information-based codebook with a codebook size of Q - [VQ.codebook]
% And generate directions for each codevector                                        - [VQ.sph]
% iter: The number of iterations for the k-means++ clustering algorithm to generate the perfect location information-based codebook
iter    = 10;                                            
VQ      = Codebook_VQ(location, BS, Q, iter);

% Save the perfect location information-based codebook
if ~isfolder('codebook')
    mkdir('codebook')
    addpath('codebook')
end
save(['codebook/VQ_codebook_H', num2str(BS.H), 'V', num2str(BS.V), '_Q', num2str(Q)], 'VQ')