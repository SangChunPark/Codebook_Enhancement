clear, clc, close all
tic
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
if ~isfolder('Figure_3/channel')
    mkdir('Figure_3/channel')
    addpath('Figure_3/channel')
end
save(['Figure_3/channel/H', num2str(BS.H), 'V', num2str(BS.V), '_sigmaL_', num2str(ch.sigma_L)], 'channel')

channel = channel./sqrt(sum(abs(channel).^2, 2));
len_ch  = length(channel);
channel2= channel;

%% Generate a DFT codebook
% Set parameters
H  = 8;                 % The number of horizontal elements 
V  = 8;                 % The number of vertical elements
N1 = H;                 % The number of horizontal elements 
N2 = V;                 % The number of vertical elements
O1 = 1;                 % The horizontal oversampling factor
O2 = 1;                 % The vertical oversampling factor
Q  = N1*N2*O1*O2;       % Codebook size

% Save the DFT codebook
if ~isfolder('Figure_3/codebook')
    mkdir('Figure_3/codebook')
    addpath('Figure_3/codebook')
end

% Generate the DFT codebook with a codebook size of Q = N1*N2*O1*O2 - [DFT.codebook]
% And generate directions for each DFT codevector                   - [DFT.sph]
% Reference: 3GPP TS 38.214 5.2.2.2 Precoding matrix indicator (PMI)
DFT = Codebook_DFT(N1, N2, O1, O2);

save(['Figure_3/codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'DFT')
initial_codebook        = DFT.codebook;
initial_codebook_sph    = DFT.sph;

%% Search parameters for codebook enhancement
% Remove codevectors which are non-unidirectivity in the initial codebook
ind_rem = (isnan(sum(initial_codebook_sph)) | isinf(sum(initial_codebook_sph)));
for i = 1 : Q
    if ~isreal(initial_codebook_sph(:, i))
        ind_rem(i) = 1;
    end
end
initial_codebook(:, ind_rem) = [];
initial_codebook_sph(:, ind_rem) = [];

% Parameter setup for kernel density estimation
alpha_cand     = 0.1 : 0.2 : 0.7;           % Smooth parameter for KDE
beta_cand      = 2.^(2 : 4);                % Initial bandwidth of KDE

margin    = 0.1;
gap       = 0.5;

KDE_setup

num_search = 0;
max_search = 5;
while num_search < max_search

    max_iteration   = 6;

    for beta = beta_cand
        for alpha = alpha_cand

            ch_ind  = randperm(len_ch, 1e3);
            channel = channel2(ch_ind, :);
            iteration = 0;

            BF_gain = zeros(1, max_iteration + 1);
            BF_gain(iteration + 1) = mean(max(abs(channel*initial_codebook).^2, [], 2));

            % Codebook Enhancement
            % Codebook Enhancement (Step 1 to Step 3)

            while iteration < max_iteration

                % STEP 1. Collect PMI using the initial codebook and convert PMI indices to UE directions

                noise = zeros(size(channel, 1), 1);

                if iteration == 0
                    codebook        = initial_codebook;
                    codebook_sph    = initial_codebook_sph;
                else
                    codebook        = enhanced.codebook;
                    codebook_sph    = enhanced.sph;
                end

                [~, PMI_list] = max(abs(channel*codebook + noise).^2, [], 2);
                x = transpose(codebook_sph(:, PMI_list));

                % STEP 2. Estimate user distribution using kernel density estimation (KDE) based on converted discrete UE directions

                bandwidth = beta*exp(-alpha*iteration);

                KDE = ksdensity(x, xi, 'Bandwidth', bandwidth);

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

                iteration = iteration + 1;

                [max_val, ~] = max(abs(channel2*enhanced.codebook).^2, [], 2);
                BF_gain(iteration + 1) = mean(max_val);
            end

            % Save beamforming gains according to the number of iterations
            if ~isfolder('Figure_3/parameter_search')
                mkdir('Figure_3/parameter_search')
                addpath('Figure_3/parameter_search')
            end

            if isfile(['Figure_3/parameter_search/BFgain_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_alpha_', num2str(alpha*10), '_beta_', num2str(beta), '.mat'])
                tmp_BF_gain = BF_gain;
                load(['Figure_3/parameter_search/BFgain_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_alpha_', num2str(alpha*10), '_beta_', num2str(beta), '.mat'])
                BF_gain = [tmp_BF_gain; BF_gain]
                save(['Figure_3/parameter_search/BFgain_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_alpha_', num2str(alpha*10), '_beta_', num2str(beta), '.mat'], 'BF_gain', 'max_iteration')
            else
                save(['Figure_3/parameter_search/BFgain_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_alpha_', num2str(alpha*10), '_beta_', num2str(beta), '.mat'], 'BF_gain', 'max_iteration')
            end

        end
    end

    num_search = num_search + 1;
end

%% Plot results of searching parameters
addpath('Figure_3/parameter_search')

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

if ~(exist('H', 'var') && exist('V', 'var') && exist('Q', 'var'))
    H = 8;                   % The number of horizontal elements
    V = 8;                   % The number of vertical elements
    Q = 64;                  % Codebook size
end

if ~(exist('alpha_cand', 'var') && exist('beta_cand', 'var'))
    alpha_cand     = 0.1 : 0.2 : 0.7;
    beta_cand      = 2.^(2 : 4);
end
len_alpha   = length(alpha_cand);
len_beta    = length(beta_cand);

figure(1)
hold on
for i = 1 : len_alpha
    plot(-1e6, -1e6, [marker{i}, 'k'], 'LineWidth', linewidth, 'MarkerSize', markersize)
end
for i = 1 : len_beta
    plot(-1e6, -1e6, linestyle{i}, 'Color', color(i, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
end

beta_cnt = 1;
for beta = beta_cand
    alpha_cnt = 1;
    for alpha = alpha_cand
        load(['Figure_3/parameter_search/BFgain_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q), '_alpha_', num2str(alpha*10), '_beta_', num2str(beta), '.mat'])
        plot(0 : max_iteration, mean(BF_gain, 1), [marker{alpha_cnt}, linestyle{beta_cnt}], 'Color', color(beta_cnt, :), 'LineWidth', linewidth, 'MarkerSize', markersize)
        alpha_cnt = alpha_cnt + 1;
    end
    beta_cnt = beta_cnt + 1;
end

legend_ = cell(1, len_alpha + len_beta);
for i = 1 : len_alpha
    legend_{i} = ['   \alpha = ', num2str(alpha_cand(i))];
end
for i = 1 : len_beta
    legend_{len_alpha + i} = ['   \beta = ', num2str(beta_cand(i))];
end
legend(legend_, 'Location' , 'best')

xlim([0 max_iteration])
% ylim([0.55 0.95])
% yticks(0.55 : 0.1 : 0.95)

grid
xlabel('The number of iterations')
ylabel('Beamforming gain')

set(gca,'FontSize', 14, 'FontName', 'Arial')
toc