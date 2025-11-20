function rxDemod = ofdmDemodulator(rxSymb, modType)
% rxSymb: subcarrier symbols as a matrix Nsc x numSymbols (as in your code)
% rxDemod: bit matrix with rows = bits_per_symbol, columns = #symbols (so flattening matches original data order)
  % column vector of received symbols
rxSymb = rxSymb(:);
switch upper(modType)
    case 'BPSK'
        bits = real(rxSymb) < 0;   % 1 when negative, 0 when positive
        rxDemod = bits;

    case 'QPSK'
        % Decision on real/imag parts
        rxScaled = rxSymb * sqrt(2);
        Ipos = real(rxScaled) > 0;
        Qpos = imag(rxScaled) > 0;
        % Recall modulator used I = 1-2*b1; Q = 1-2*b2
        % So b1 = (I<0), b2 = (Q<0)
        b1 = ~Ipos;  % 1 when I <=0
        b2 = ~Qpos;
        % Compose in same order used in modulator: [b1 b2] with left-msb
        bits = [b1 b2];   % Nx2
        rxDemod = bits.'; % 2 x N (matches reshape->(:) flatten order)
        rxDemod = rxDemod(:);

    case 'ASK'
        % original levels used by the modulator
    levels = [-3 -1 1 3];    % 1x4
    % apply the same normalization the modulator used
    normFactor = sqrt(mean(abs(levels).^2));
    levels_norm = levels / normFactor;   % 1x4, matches tx scaling

    % ensure column
    rxSymb = rxSymb(:);
    % distances: Nx4
    dist = abs(rxSymb - levels_norm);
    [~, idx] = min(dist, [], 2);  % idx in 1..4

    % convert back to bits (binary, left-msb) - matches your modulator
    bits = de2bi(idx-1, 2, 'left-msb'); % Nx2
    rxDemod = bits.'; 
    rxDemod = rxDemod(:);

    case 'FSK'
        M = 4;
        k = 0:(M-1);
        ref = exp(1j*2*pi*k/M);   % 1xM
        dist = abs(rxSymb - ref);  % NxM
        [~, idx] = min(dist, [], 2);
        bits = de2bi(idx-1, log2(M), 'left-msb'); % Nx2
        rxDemod = bits.'; rxDemod = rxDemod(:);

    case '8-PSK'
        M = 8;
        % pskdemod expects symbols and returns integer symbol indices
        idx = pskdemod(rxSymb, M, pi/M);  % column of 0..M-1
        bits = de2bi(idx, log2(M), 'left-msb'); % Nx3
        rxDemod = bits.'; rxDemod = rxDemod(:);

    case {'16QAM'}
        % Use same options as modulator
        idx = qamdemod(rxSymb, 16, 'gray', 'UnitAveragePower', true);
        bits = de2bi(idx, 4, 'left-msb');  % Nx4
        rxDemod = bits.'; rxDemod = rxDemod(:);

    otherwise
        error('Unsupported modulation type: %s', modType);
end
end
