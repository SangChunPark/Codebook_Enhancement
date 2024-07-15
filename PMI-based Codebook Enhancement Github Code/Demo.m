clear, clc, close all

addpath('function')

if ~isfolder('Demo')
    mkdir('Demo')
    addpath('Demo')
end

default_path = 'Demo/';

%% Environment
tic
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
n_ch            = 4e5;
load('location_xs.mat')                 % location (x_UE, y_UE) [m], dimensions:(# of UE) x 2
location(n_ch + 1 : end, :) = [];
location(:, 3)  = UE.h;                 % location (x_UE, y_UE, h_UE) [m], dimensions: (# of UE) x 3
len_channel     = size(location, 1);
location        = transpose(location);  % Dimensions: 3 x (# of UE)

% Generate channels using Saleh-Valenzuela (SV) channel model
ch.L           = 20;               % The number of multi-paths
ch.sigma_alpha = 1;                % Std. dev. of complex gain for channels
ch.sigma_L     = 10;                % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
ch.radiation_pattern = false;      % Radiation power pattern: ture - The radiation power pattern based on 3GPP TR 38.901, false - Isotropic power pattern

t = toc;
disp(['Set the environment: ', num2str(round(t, 3)), ' sec.'])
%% Channel Generation
tic
channel = zeros(len_channel, BS.H*BS.V);
for i = 1 : size(location, 2)
    UE.loc = location(:, i);
    channel(i, :) = SV_channel(BS, UE, ch);
end

% Save the channel
if ~isfolder([default_path, 'channel'])
    mkdir([default_path, 'channel'])
    addpath([default_path, 'channel'])
end
save([default_path, 'channel/H', num2str(BS.H), 'V', num2str(BS.V), '_sigmaL_', num2str(ch.sigma_L)], 'channel')

t = toc;
disp(['Generate ', num2str(n_ch), ' channels with std. of Laplacian distribution ', num2str(ch.sigma_L), ' [deg] : ', num2str(round(t, 3)), ' sec.'])
%% VQ Codebook Generation
tic
Q       = 64;              % Codebook size

% Remove duplicate location coordinates to reduce computational load in function of "Codebook_VQ".
if size(location, 2) ~= 3
    location        = transpose(location);
end
location        = round(location, 3);

% Generate the perfect location information-based codebook with a codebook size of Q - [VQ.codebook]
% And generate directions for each codevector                                        - [VQ.sph]
% iter: The number of iterations for the k-means++ clustering algorithm to generate the perfect location information-based codebook
iter    = 10;
VQ      = Codebook_VQ(location, BS, Q, iter);

% Save the perfect location information-based codebook
if ~isfolder([default_path, 'codebook'])
    mkdir([default_path, 'codebook'])
    addpath([default_path, 'codebook'])
end
save([default_path, 'codebook/VQ_codebook_H', num2str(BS.H), 'V', num2str(BS.V), '_Q', num2str(Q)], 'VQ')

t = toc;
disp(['Generate a VQ codebook with a codebook size of ', num2str(Q), ': ', num2str(round(t, 3)), ' sec.'])
%% DFT Codebook Generation
tic
% Set parameters
H  = BS.H;                 % The number of horizontal elements
V  = BS.V;                 % The number of vertical elements
N1 = H;                    % The number of horizontal elements
N2 = V;                    % The number of vertical elements
O1 = 1;                    % The horizontal oversampling factor
O2 = 1;                    % The vertical oversampling factor
Q  = N1*N2*O1*O2;          % Codebook size

% Generate the DFT codebook with a codebook size of Q = N1*N2*O1*O2 - [DFT.codebook]
% And generate directions for each DFT codevector                   - [DFT.sph]
% Reference: 3GPP TS 38.214 5.2.2.2 Precoding matrix indicator (PMI)
DFT = Codebook_DFT(N1, N2, O1, O2);

% Save the DFT codebook
save([default_path, 'codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'DFT')

t = toc;
disp(['Generate a DFT codebook with a codebook size of ', num2str(Q), ': ', num2str(round(t, 3)), ' sec.'])
%% Codebook Enhancement from the DFT Codebook
tic

% Initial codebook dimensions: (# of BS antenna elements) X (codebook size)
% Initial codebook directions (theta, phi), dimensions: 2 X (codebook size)
initial_codebook     = DFT.codebook;
initial_codebook_sph = DFT.sph;

% Remove codevectors which are non-unidirectivity in the initial codebook
ind_rem = (isnan(sum(initial_codebook_sph)) | isinf(sum(initial_codebook_sph)));
for i = 1 : Q
    if ~isreal(initial_codebook_sph(:, i))
        ind_rem(i) = 1;
    end
end
initial_codebook(:, ind_rem)        = [];
initial_codebook_sph(:, ind_rem)    = [];

% Parameter setup for kernel density estimation
alpha     = 0.1;        % Smooth parameter for KDE
beta      = 4;          % Initial bandwidth of KDE

margin    = 0.1;
gap       = 0.5;

KDE_setup

% Codebook Enhancement (Step 1 to Step 3)

max_iteration   = 6;
iteration       = 0;

while iteration < max_iteration

    ch_ind  = randperm(length(channel), 1e3);
    channel2 = channel(ch_ind, :);
    % STEP 1. Collect PMI using the initial codebook and convert PMI indices to UE directions

    noise = zeros(size(channel2, 1), 1);

    if iteration == 0
        codebook        = initial_codebook;
        codebook_sph    = initial_codebook_sph;
    else
        codebook        = enhanced.codebook;
        codebook_sph    = enhanced.sph;
    end

    [~, PMI_list] = max(abs(channel2*codebook + noise).^2, [], 2);
    x = transpose(codebook_sph(:, PMI_list));

    % Plot1

    % STEP 2. Estimate user distribution using kernel density estimation (KDE) based on converted discrete UE directions

    bandwidth = beta*exp(-alpha*iteration);

    KDE = ksdensity(x, xi, 'Bandwidth', bandwidth);

    % Plot2

    % STEP 3. Determine k centroids using the k-means++ clustering based on estimated user distribution

    KDE2 = zeros(length(xj), 1);
    for i = 1 : length(ic)
        KDE2(ic(i)) = KDE2(ic(i)) + KDE(i);
    end

    scale = 1e5;
    KDE_scaled = round(KDE2*scale);
    sample = xj(KDE_scaled ~= 0, :);
    KDE_scaled = KDE_scaled(KDE_scaled ~= 0);

    iter = 5;              % The number of iterations for the k-means clustering algorithm
    enhanced.sph = kmeans_plus_clustering(sample, KDE_scaled, Q, iter);
    enhanced.codebook = zeros(H*V, Q);
    for i = 1 : Q
        enhanced.codebook(:, i) = ARV_UPA(enhanced.sph(1, i), enhanced.sph(2, i), H, V)./sqrt(H*V);
    end

    % Plot3

    iteration = iteration + 1;

    [max_val, ~] = max(abs(channel*enhanced.codebook).^2, [], 2);
    BF_gain(iteration + 1) = mean(max_val);
end


% Save the enhanced codebook
save([default_path, 'codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'enhanced')

t = toc;
disp(['Generate a enhanced codebook with a codebook size of ', num2str(Q), ' from the DFT codebook: ', num2str(round(t, 3)), ' sec.'])
%% Performance Evaluation Sum-Rate [bps/Hz]
tic
codebook_cand = {'VQ', 'enhanced', 'DFT'};
channel = channel./sqrt(sum(abs(channel).^2, 2));

for i = 1 : length(codebook_cand)

    load([default_path, 'codebook/', codebook_cand{i}, '_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
    switch codebook_cand{i}
        case 'enhanced'
            codebook = enhanced.codebook;
        case 'VQ'
            codebook = VQ.codebook;
        case 'DFT'
            codebook = DFT.codebook;
    end

    % Rank adaptive sum-rate
    noise_power_dB    = -10;
    noise_power       = 10^(noise_power_dB/10);
    len_channel       = size(channel, 1);
    sum_rate          = zeros(H*V, 2);

    max_iteration = 1e5;
    sum_rate_data = zeros(2, max_iteration);
    warning ('off','all')

    rank_adaptive_sum_rate

    if ~isfolder([default_path, 'sum_rate'])
        mkdir([default_path, 'sum_rate'])
        addpath([default_path, 'sum_rate'])
    end
    save([default_path, 'sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(ch.sigma_L), '_', codebook_cand{i}], 'sum_rate')


    sum_rate(sum_rate(:, 2) < max_iteration/1000, 2) = 0;

    plot(sum_rate(:, 1)./sum_rate(:, 2), 'LineWidth', 2)

    disp(['Evaluate sum-rate using the ',  codebook_cand{i}, ' codebook, the average sum-rate = ', num2str(sum(sum_rate(:, 1))./sum(sum_rate(:, 2))), ' bps/Hz'])

    hold on
end
t = toc;
disp(['Evaluate sum-rate performance: ', num2str(round(t, 3)), ' sec.'])

xlim([0 15])
ylim([0 20])

grid
xlabel('Rank')
ylabel('Sum-rate [bps/Hz]')

legend('   Perfect location information', '   Enhanced DFT codebook', '   Conventional DFT codebook', '    ', '    ', '    ', 'Location' , 'NorthWest')
legend boxoff

set(gca,'FontSize', 14, 'FontName', 'Arial')