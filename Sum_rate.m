clear, clc, close all

% Environment Setup

addpath('function', 'channel', 'codebook')

H = 8;                  % The number of horizontal elements 
V = 8;                  % The number of vertical elements
Q = 256;                 % Codebook size

% Load channels, channel dimensions: (# of UE) X (# of BS antenna elements)
sigma_L = 10;
load(['channel/H', num2str(H), 'V', num2str(V), '_sigmaL_', num2str(sigma_L)]);
channel = channel./sqrt(sum(abs(channel).^2, 2));

% Load the codebook
type_codebook = 'DFT';
load([type_codebook, '_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)])
switch type_codebook
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

if ~isfolder('sum_rate')
    mkdir('sum_rate')
    addpath('sum_rate')
end
save(['sum_rate/H', num2str(H), 'V', num2str(V), '_Q', num2str(Qc), '_sigmaL_', num2str(sigma_L), '_', type_codebook], 'sum_rate')