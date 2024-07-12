addpath('function')

% Set parameters
H  = 8;                 % The number of horizontal elements 
V  = 8;                 % The number of vertical elements
N1 = H;                 % The number of horizontal elements 
N2 = V;                 % The number of vertical elements
O1 = 1;                 % The horizontal oversampling factor
O2 = 1;                 % The vertical oversampling factor
Q  = N1*N2*O1*O2;       % Codebook size

% Generate the DFT codebook with a codebook size of Q = N1*N2*O1*O2 - [DFT.codebook]
% And generate directions for each DFT codevector                   - [DFT.sph]
% Reference: 3GPP TS 38.214 5.2.2.2 Precoding matrix indicator (PMI)
DFT = Codebook_DFT(N1, N2, O1, O2);

% Save the DFT codebook
if ~isfolder('codebook')
    mkdir('codebook')
    addpath('codebook')
end
save(['codebook/DFT_codebook_H', num2str(H), 'V', num2str(V), '_Q', num2str(Q)], 'DFT')