clear, clc, close all

% Environment Setup

addpath('function', 'channel', 'codebook')

H = 8;                  % The number of horizontal elements 
V = 8;                  % The number of vertical elements
Q = 64;                 % Codebook size

% Load channels, channel dimensions: (# of UE) X (# of BS antenna elements)
sigma_L = 0;
load(['channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)]);
channel = channel./sqrt(sum(abs(channel).^2, 2));

% Initial codebook dimensions: (# of BS antenna elements) X (codebook size)
% Initial codebook directions (theta, phi), dimensions: 2 X (codebook size)
load(['DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
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

max_iteration   = 5;
iteration       = 0;
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
if ~isfolder('codebook')
    mkdir('codebook')
    addpath('codebook')
end
save(['codebook/enhanced_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'enhanced')