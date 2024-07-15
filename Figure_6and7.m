clear, clc, close all

addpath('function', 'channel', 'codebook')

%% Generate channels
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
location(4e4 + 1 : end, :) = [];
location(:, 3)  = UE.h;                 % location (x_UE, y_UE, h_UE) [m], dimensions: (# of UE) x 3
len_channel     = size(location, 1);
location        = transpose(location);  % Dimensions: 3 x (# of UE)

% Generate channels using Saleh-Valenzuela (SV) channel model
ch.L           = 20;               % The number of multi-paths
ch.sigma_alpha = 1;                % Std. dev. of complex gain for channels
ch.radiation_pattern = false;      % Radiation power pattern: ture - The radiation power pattern based on 3GPP TR 38.901, false - Isotropic power pattern

% Save the channel
if ~isfolder('Figure_6&7/channel')
    mkdir('Figure_6&7/channel')
    addpath('Figure_6&7/channel')
end

for sigma_L = 0 : 5 : 10           % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
    ch.sigma_L = sigma_L;
    channel = zeros(len_channel, BS.H*BS.V);
    for i = 1 : size(location, 2)
        UE.loc = location(:, i);
        channel(i, :) = SV_channel(BS, UE, ch);
    end

    save(['Figure_6&7/channel/H', num2str(BS.H), 'V', num2str(BS.V), '_sigmaL_', num2str(sigma_L)], 'channel')
end

channel = channel./sqrt(sum(abs(channel).^2, 2));
len_ch  = length(channel);
channel2= channel;

%% Generate DFT codebooks
% Set parameters
H  = 8;                 % The number of horizontal elements
V  = 8;                 % The number of vertical elements
N1 = H;                 % The number of horizontal elements
N2 = V;                 % The number of vertical elements

% Save the DFT codebook
if ~isfolder('Figure_6&7/codebook')
    mkdir('Figure_6&7/codebook')
    addpath('Figure_6&7/codebook')
end

for O_tmp = [1/2, 1, 2, 4]
    O1 = O_tmp;
    O2 = O_tmp;
    Q  = N1*N2*O1*O2;       % Codebook size
    % Generate the DFT codebook with a codebook size of Q = N1*N2*O1*O2 - [DFT.codebook]
    % And generate directions for each DFT codevector                   - [DFT.sph]
    % Reference: 3GPP TS 38.214 5.2.2.2 Precoding matrix indicator (PMI)
    DFT = Codebook_DFT(N1, N2, O1, O2);

    save(['Figure_6&7/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'DFT')
    initial_codebook        = DFT.codebook;
    initial_codebook_sph    = DFT.sph;
end

%% Generate VQ codebooks

% Load UE locations
load('location_xs.mat')    % location (x_UE, y_UE) [m], dimensions:(# of UE) x 2
location(5e4 + 1 : end, :) = [];
location(:, 3)  = UE.h;     % location (x_UE, y_UE, h_UE) [m], dimensions: (# of UE) x 3


% Remove duplicate location coordinates to reduce computational load in function of "Codebook_VQ".
location        = round(location, 3);

iter    = 10;
for Q = 2.^(4 : 2 : 10)
    % Generate the perfect location information-based codebook with a codebook size of Q - [VQ.codebook]
    % And generate directions for each codevector                                        - [VQ.sph]
    % iter: The number of iterations for the k-means++ clustering algorithm to generate the perfect location information-based codebook
    VQ      = Codebook_VQ(location, BS, Q, iter);

    save(['Figure_6&7/codebook/VQ_codebook_H', num2str(BS.H), 'V', num2str(BS.V), '_Q', num2str(Q)], 'VQ')
end

%% Search parameters for codebook enhancement

% Parameter setup for kernel density estimation
alpha     = 0.1;        % Smooth parameter for KDE
beta      = 4;          % Initial bandwidth of KDE

margin    = 0.1;
gap       = 0.5;

KDE_setup

for sigma_L = [0, 10]   % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
    load(['Figure_6&7/channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)])

    for Q = 2.^(4 : 2 : 10)
        load(['Figure_6&7/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
        initial_codebook        = DFT.codebook;
        initial_codebook_sph    = DFT.sph;

        % Remove codevectors which are non-unidirectivity in the initial codebook
        ind_rem = (isnan(sum(initial_codebook_sph)) | isinf(sum(initial_codebook_sph)));
        for i = 1 : Q
            if ~isreal(initial_codebook_sph(:, i))
                ind_rem(i) = 1;
            end
        end
        initial_codebook(:, ind_rem) = [];
        initial_codebook_sph(:, ind_rem) = [];

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
        end
        save(['Figure_6&7/codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L)], 'enhanced')
    end
end

sigma_L = 5;
load(['Figure_6&7/channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)])
Q = 64;
load(['Figure_6&7/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
initial_codebook        = DFT.codebook;
initial_codebook_sph    = DFT.sph;

% Remove codevectors which are non-unidirectivity in the initial codebook
ind_rem = (isnan(sum(initial_codebook_sph)) | isinf(sum(initial_codebook_sph)));
for i = 1 : Q
    if ~isreal(initial_codebook_sph(:, i))
        ind_rem(i) = 1;
    end
end
initial_codebook(:, ind_rem) = [];
initial_codebook_sph(:, ind_rem) = [];

margin    = 0.1;
gap       = 0.5;

KDE_setup

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
end
save(['Figure_6&7/codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L)], 'enhanced')

%% Evaluate sum-rate performance
warning ('off','all')

noise_power_dB    = -10;
noise_power       = 10^(noise_power_dB/10);
max_iteration     = 1e5;

if ~isfolder('Figure_6&7/sum_rate')
    mkdir('Figure_6&7/sum_rate')
    addpath('Figure_6&7/sum_rate')
end

sigma_cnt = 1;
for sigma_L = [0, 10]   % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
    load(['Figure_6&7/channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)])
    channel = channel./sqrt(sum(abs(channel).^2, 2));
    len_channel       = size(channel, 1);
    Q_cnt = 1;
    for Q = 2.^(4 : 2 : 10)
        load(['Figure_6&7/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
        load(['Figure_6&7/codebook/VQ_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
        load(['Figure_6&7/codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L)])

        % Rank adaptive sum-rate
        codebook          = DFT.codebook;
        sum_rate          = zeros(H*V, 2);
        sum_rate_data     = zeros(2, max_iteration);

        rank_adaptive_sum_rate
        save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_DFT'], 'sum_rate', 'sum_rate_data')

        % Rank adaptive sum-rate
        codebook          = VQ.codebook;
        sum_rate          = zeros(H*V, 2);
        sum_rate_data     = zeros(2, max_iteration);

        rank_adaptive_sum_rate
        save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_VQ'], 'sum_rate', 'sum_rate_data')

        % Rank adaptive sum-rate
        codebook          = enhanced.codebook;
        sum_rate          = zeros(H*V, 2);
        sum_rate_data     = zeros(2, max_iteration);

        rank_adaptive_sum_rate
        save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_enhanced'], 'sum_rate', 'sum_rate_data')

        Q_cnt = Q_cnt + 1;
    end
    sigma_cnt = sigma_cnt + 1;
end

sigma_L = 5;   % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
load(['Figure_6&7/channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)])
channel = channel./sqrt(sum(abs(channel).^2, 2));
len_channel       = size(channel, 1);

Q = 64;
load(['Figure_6&7/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
load(['Figure_6&7/codebook/VQ_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
load(['Figure_6&7/codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L)])

% Rank adaptive sum-rate
codebook          = DFT.codebook;
sum_rate          = zeros(H*V, 2);
sum_rate_data     = zeros(2, max_iteration);

rank_adaptive_sum_rate
save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_DFT'], 'sum_rate', 'sum_rate_data')

% Rank adaptive sum-rate
codebook          = VQ.codebook;
sum_rate          = zeros(H*V, 2);
sum_rate_data     = zeros(2, max_iteration);

rank_adaptive_sum_rate
save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_VQ'], 'sum_rate', 'sum_rate_data')

% Rank adaptive sum-rate
codebook          = enhanced.codebook;
sum_rate          = zeros(H*V, 2);
sum_rate_data     = zeros(2, max_iteration);

rank_adaptive_sum_rate
save(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_enhanced'], 'sum_rate', 'sum_rate_data')

%% Plot performance evaluations
clear, clc, close all

H = 8;
V = 8;

color(1, :) = [0, 0.4470, 0.7410];
color(2, :) = [0.8500, 0.3250, 0.0980];
color(3, :) = [0.9290, 0.6940, 0.1250];
color(4, :) = [0.4660, 0.6740, 0.1880];
color(5, :) = [0.4940, 0.1840, 0.5560];
color(6, :) = [0.3010 0.7450 0.9330];

linewidth   = 2;
markersize  = 10;

marker      = {'o', 's', '^', 'x', '>', 'd'};
linestyle   = {'-', '--', '-.', ':'};

mean_SR_DFT         = zeros(3, 1);
mean_SR_VQ          = zeros(3, 1);
mean_SR_enhanced    = zeros(3, 1);

Q = 64;
sigma_cnt = 1;
for sigma_L = [0, 5, 10]   % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]

    load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_DFT'])
    mean_SR_DFT(sigma_cnt, 1) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

    load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_VQ'])
    mean_SR_VQ(sigma_cnt, 1) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

    load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_enhanced'])
    mean_SR_enhanced(sigma_cnt, 1) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

    sigma_cnt = sigma_cnt + 1;
end

figure(1)
hold on

for i = [2, 1, 3]
    plot(-1e6, -1e6, linestyle{i}, 'Color', color(i, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
end

plot(0 : 5 : 10, mean_SR_VQ, linestyle{2}, 'Color', color(2, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
plot(0 : 5 : 10, mean_SR_enhanced, linestyle{1}, 'Color', color(1, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
plot(0 : 5 : 10, mean_SR_DFT, linestyle{3}, 'Color', color(3, :), 'LineWidth', linewidth, 'MarkerSize', markersize)

xlim([-2.5 12.5])
xticks(0 : 5 : 10)
ylim([0 20])
yticks(0 : 5 : 20)

grid
xlabel('Angle spread [deg]')
ylabel('Sum-rate [bps/Hz]')

legend('   Perfect location information', '   Enhanced DFT codebook', '   Conventional DFT codebook', 'Location' , 'NorthWest')
legend boxoff

set(gca,'FontSize', 14, 'FontName', 'Arial')


mean_SR_DFT         = zeros(2, 4        );
mean_SR_VQ          = zeros(2, 4);
mean_SR_enhanced    = zeros(2, 4);

sigma_cnt = 1;
for sigma_L = [0, 10]   % Std. dev. of Laplacian distribution for azimuth/zenith angles of multi-paths [deg]
    Q_cnt = 1;
    for Q = 2.^(4 : 2 : 10)
        load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_DFT'])
        mean_SR_DFT(sigma_cnt, Q_cnt) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

        load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_VQ'])
        mean_SR_VQ(sigma_cnt, Q_cnt) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

        load(['Figure_6&7/sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_sigmaL_', num2str(sigma_L), '_enhanced'])
        mean_SR_enhanced(sigma_cnt, Q_cnt) = sum(sum_rate(:, 1))./sum(sum_rate(:, 2));

        Q_cnt = Q_cnt + 1;
    end
    sigma_cnt = sigma_cnt + 1;
end

figure(2)
hold on

for i = [2, 1, 3]
    plot(-1e6, -1e6, linestyle{i}, 'Color', color(i, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
end

for i = [1, 2]
    plot(-1e6, -1e6, [marker{i}, 'k'], 'LineWidth', linewidth, 'MarkerSize', markersize)
end

for i = [1, 2]
plot(4 : 2 : 10, mean_SR_VQ(i, :), [linestyle{2}, marker{i}], 'Color', color(2, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
plot(4 : 2 : 10, mean_SR_enhanced(i, :), [linestyle{1}, marker{i}], 'Color', color(1, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
plot(4 : 2 : 10, mean_SR_DFT(i, :), [linestyle{3}, marker{i}], 'Color', color(3, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
end

xlim([3 11])
xticks(4 : 2 : 10)
ylim([0 25])
yticks(0 : 5 : 25)

grid
xlabel('The number of feedback bits {\itB}')
ylabel('Sum-rate [bps/Hz]')

legend('   Perfect location information', '   Enhanced DFT codebook', '   Conventional DFT codebook', 'Location' , 'NorthWest')
legend boxoff

set(gca,'FontSize', 14, 'FontName', 'Arial')