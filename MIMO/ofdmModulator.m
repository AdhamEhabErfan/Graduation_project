function txMod = ofdmModulator(data, modType)
% data : column vector of bits (0/1)
% txMod: column vector of complex symbols
data = data(:); % ensure column

switch upper(modType)
    case 'BPSK'
        % 1 bit/symbol: 0 -> +1, 1 -> -1
        txMod = 1 - 2*data;       % +1 / -1
        txMod = txMod(:);

    case 'QPSK'
        bitsPerSym = 2;
        B = reshape(data, bitsPerSym, []).'; % rows: symbols, cols: bits (MSB left)
        % Gray mapping: 00 ->  1+1j, 01 ->  1-1j, 11 -> -1-1j, 10 -> -1+1j
        I = 1 - 2*B(:,1); % +1 / -1
        Q = 1 - 2*B(:,2); % +1 / -1
        % Reorder Q based on Gray mapping (flip sign for second bit)
        % To match demod code below we use (I) + 1j*(Q) and then normalize
        tx = (I + 1j*Q) ./ sqrt(2);
        txMod = tx(:);

    case 'ASK'  % 4-PAM (2 bits/symbol)
        bitsPerSym = 2;
        B = reshape(data, bitsPerSym, []).';            % Nx2
        idx = bi2de(B, 'left-msb');                     % 0..3
        levels = [-3 -1 1 3];                           % 4-PAM levels
        tx = levels(idx+1).';
        % normalize to unit average power
        tx = tx / sqrt(mean(abs(tx).^2));
        txMod = tx(:);

    case 'FSK'  % M-FSK, orthogonal tones represented as complex phasors (symbolic)
        M = 4;
        bitsPerSym = log2(M);
        B = reshape(data, bitsPerSym, []).';
        idx = bi2de(B, 'left-msb'); % 0..M-1
        % Use discrete orthogonal phasors (equal-energy)
        k = 0:(M-1);
        ref = exp(1j*2*pi*k/M);            % 1xM
        tx = ref(idx+1).';                  % column vector
        txMod = tx(:);

    case '8-PSK'  % M-PSK, use built-in to avoid mistakes
        M = 8;
        bitsPerSym = log2(M);
        B = reshape(data, bitsPerSym, []).';
        idx = bi2de(B, 'left-msb');           % 0..M-1
        % Use pskmod with phase offset pi/M to get Gray-like mapping consistent with demod
        tx = pskmod(idx, M, pi/M);
        txMod = tx(:);

    case {'16QAM'}
        bitsPerSym = 4;
        B = reshape(data, bitsPerSym, []).';
        idx = bi2de(B, 'left-msb');           % 0..15
        % Use unit average power so scaling consistent with demod
        tx = qammod(idx, 16, 'gray', 'UnitAveragePower', true);
        txMod = tx(:);

    otherwise
        error('Unsupported modulation type: %s', modType);
end
end
