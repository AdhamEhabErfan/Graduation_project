rng(1); % Consistent randomness
%====================================================================
% Assuming we have Bandwidth of 25 MHz
% and subcarrier spacing (delta f) of 30 kHz.

   % Number of subcarriers= 833.33. It is rounded to 840 (multiple of 12)
   %In 4G/5G: 1 RB (Resource Block) = 12 consecutive subcarriers per slot in time
   % 840(70 RBs x 12)

      %FFT size is set to be the nearest power of two(2^x) next to Number of subcarriers
      %let FFT size (NFFT) be 2^10 = 1024
      %sampling frequenccy fs= NFFT * delta f = 30.72 MHz
      %cyclic prefix length = cyclic prefix time/ sampling time
      % = 2.38us * fs =~ 73

DeltaF = 30e3;   % Subcarrier spacing
Nsc    = 840;    % Number of subcarriers
NFFT   = 1024;   % FFT size
cpLen  = 73;     % Cyclic prefix length
fs     = NFFT * DeltaF;

N_tx = 2;        % Transmit antennas
N_rx = 2;        % Receive antennas

totalBits = 1e5;
SNRdB = 0:5:20; 
modSchemes = {'BPSK'}; % Extend as needed

mimoMode = 'multiplexing'; % 'multiplexing' or 'diversity'

figure; hold on; grid on;

for mIdx = 1:length(modSchemes)
    modType = modSchemes{mIdx};

    switch modType
        case 'BPSK'
            bitsPerSymbol = 1;
        case 'QPSK'
            bitsPerSymbol = 2;
        case 'ASK'
            bitsPerSymbol = 2;
        case 'FSK'
            bitsPerSymbol = 2;
        case '8-PSK'
            bitsPerSymbol = 3;
        case '16QAM'
            bitsPerSymbol = 4;
    end

    
    if strcmp(mimoMode,'multiplexing')
    % Multiplexing: throughput divided by N_tx
    numOFDMSymbols = floor(totalBits/(Nsc * bitsPerSymbol * N_tx)); % number of OFDM symbols per transmitting antenna per sub-carrier
    numSymbolsPerStream = Nsc * numOFDMSymbols;
    dataBitsPerTx = numSymbolsPerStream * bitsPerSymbol;
        dataLen = dataBitsPerTx * N_tx;
        data = randi([0 1], dataLen, 1); % Total Random bits for all antennas
        data_reshaped = reshape(data, N_tx, dataBitsPerTx);  % Each row corresponds to one Tx antenna
    else
    % Diversity: full bitstream per antenna
    numOFDMSymbols = floor(totalBits/(Nsc * bitsPerSymbol));
    numSymbolsPerStream = Nsc * numOFDMSymbols;
    dataBitsPerTx = numSymbolsPerStream * bitsPerSymbol;
        % Generate only one stream worth of data
        data = randi([0 1], dataBitsPerTx, 1);  % Generate single stream
        data_reshaped = repmat(data', N_tx, 1); % Duplicate for all antennas (rows)
    end

    % Preallocate BER vector
    BER = zeros(length(SNRdB),1);

    % Loop over SNR values
    for idx = 1:length(SNRdB)
        nErrors_total = 0;
        bits_processed = 0;
        txSymb_grid = zeros(Nsc, numOFDMSymbols, N_tx);

        % Modulate bits for each transmit antenna (or duplicate stream for diversity)
        for working_tx_antenna = 1:N_tx
            bits_this_ant = data_reshaped(working_tx_antenna, :)';
            txMod = ofdmModulator(bits_this_ant, modType); % Convert bits to symbols
            txModParallel = reshape(txMod, Nsc, numOFDMSymbols); % Arrange symbols into OFDM grid
            txSymb_grid(:, :, working_tx_antenna) = txModParallel;
        end

        % Loop over each OFDM symbol
        for symb_idx = 1:numOFDMSymbols
            % Generate MIMO channel per subcarrier in frequency domain
            H_freq = (randn(N_rx, N_tx, Nsc) + 1i * randn(N_rx, N_tx, Nsc))*(1/sqrt(2));
            rx_symbols_allSCs = zeros(N_rx, Nsc); % Preallocate received symbols (no noise yet)

            % Transmit through channel (no noise yet)
            for sc = 1:Nsc
                H_sc = H_freq(:, :, sc);                 % N_rx x N_tx
                tx_symbols_SC = squeeze(txSymb_grid(sc, symb_idx, :)); % N_tx x 1
                rx_symbols_allSCs(:, sc) = H_sc * tx_symbols_SC; % received vector (N_rx x 1)
            end

            SNRlin = 10^(SNRdB(idx)/10); % SNRlin = Es/No where Es=1 (normalized)
            noise_power = 1 / (SNRlin);
            noise = sqrt(noise_power/2) * (randn(size(rx_symbols_allSCs)) + ...
                                         1i * randn(size(rx_symbols_allSCs)));
            rx_symbols_noisy_withFading = rx_symbols_allSCs + noise;

            % Prepare containers for equalized/combined symbols
            if strcmp(mimoMode,'multiplexing')
                % For multiplexing we need N_tx estimates per subcarrier
                rxSymb_equalized = zeros(N_tx, Nsc);
                % MMSE equalization per subcarrier -> estimate each TX stream
                for sc = 1:Nsc
                    H_sc = H_freq(:, :, sc); % N_rx x N_tx
                    % MMSE equalizer: W = (H' H + N0 I)^(-1) H'
                    mmse = inv((H_sc' * H_sc + noise_power * eye(N_tx))) * H_sc'; % N_tx x N_rx
                    rxSymb_equalized(:, sc) = mmse * rx_symbols_noisy_withFading(:, sc); % N_tx x 1
                end

                % Demodulate and compute bit errors for each tx stream
                for tx_idx = 1:N_tx
                    rxSymb = rxSymb_equalized(tx_idx, :);  % 1 x Nsc for this TX stream
                    rxDemod = ofdmDemodulator(rxSymb, modType);
                    rxBits = rxDemod(:);

                    % region of bits in original data that corresponds to the current symbol(iteration)
                    offset = (symb_idx-1)*Nsc*bitsPerSymbol + 1;
                    stop = symb_idx*Nsc*bitsPerSymbol;

                    OriginalDataSlices = data_reshaped(tx_idx, offset:stop);

                    [nErrors, ~] = biterr(OriginalDataSlices(:), rxBits);
                    nErrors_total = nErrors_total + nErrors;
                    bits_processed = bits_processed + length(OriginalDataSlices);
                end

            else % Diversity mode: perform MRC combining to get single stream estimate
                combined_symbols = zeros(1, Nsc); % 1 x Nsc combined result

                for sc = 1:Nsc
                    H_sc = H_freq(:, :, sc);             % N_rx x N_tx
                    % Effective channel for replicated transmissions: sum across columns
                    h_eff = H_sc * ones(N_tx, 1);       % N_rx x 1 (each rx antenna sees sum of tx contributions)
                    y_sc = rx_symbols_noisy_withFading(:, sc); % N_rx x 1

                    % If h_eff is all zeros (very unlikely), avoid division by zero:
                    denom = sum(abs(h_eff).^2);
                    if denom == 0
                        combined_symbols(sc) = 0;
                    else
                        % MRC combining: (h_eff^H * y) / (h_eff^H * h_eff)
                        combined_symbols(sc) = (conj(h_eff).' * y_sc) / denom;
                    end
                end

                % Demodulate combined stream once
                rxDemod = ofdmDemodulator(combined_symbols, modType);
                rxBits = rxDemod(:);

                % Compare with original single stream (stored in first row of data_reshaped)
                offset = (symb_idx-1)*Nsc*bitsPerSymbol + 1;
                stop = symb_idx*Nsc*bitsPerSymbol;
                OriginalDataSlices = data_reshaped(1, offset:stop);

                [nErrors, ~] = biterr(OriginalDataSlices(:), rxBits);
                nErrors_total = nErrors_total + nErrors;
                bits_processed = bits_processed + length(OriginalDataSlices);
            end
        end       
        BER(idx) = nErrors_total / bits_processed;
        fprintf('SNR = %d dB for %s... BER = %e\n', SNRdB(idx), modType, BER(idx));
    end

    % ========================= GRAPH PLOTTING =========================
    colors = lines(2); % one color for each MIMO Mode
    if strcmp(mimoMode,'multiplexing')
        colorIdx = 1;
    else
        colorIdx = 2;
    end

    semilogy(SNRdB, BER, 'LineWidth', 2, ...
             'Color', colors(colorIdx,:), ...
             'DisplayName', [modType ' ' num2str(N_tx) 'x' num2str(N_rx) ' MIMO (' mimoMode ')']);
    % ==================================================================
end


set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title([num2str(N_tx) 'x' num2str(N_rx) ' MIMO - ' upper(mimoMode)]);
legend('show', 'Location', 'southwest');
ylim([1e-4, 1]);
grid on;