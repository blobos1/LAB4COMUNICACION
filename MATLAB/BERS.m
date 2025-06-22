% Parametros
EbN0_dB = 1:11;                       % Rango Eb/N0 
num_bits = 1e5;                       
iter_max = 10000000;                     
target_err = 100;                    

BER_BPSK = zeros(1, length(EbN0_dB));
BER_QPSK = zeros(1, length(EbN0_dB));
BER_8PSK = zeros(1, length(EbN0_dB));


% BPSK

M = 2; k = log2(M);
const_BPSK = [-1 1];

for i = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(i)/10);
    noise_var = 1 / (2 * EbN0);  

    err = 0; iter = 0;
    while err < target_err && iter < iter_max
        bits = randi([0 1], 1, num_bits);
        s = 2*bits - 1;                        
        noise = sqrt(noise_var)*randn(1, num_bits);
        r = s + noise;
        bits_hat = r > 0;
        err = err + sum(bits ~= bits_hat);
        iter = iter + 1;
    end
    BER_BPSK(i) = err / (num_bits * iter);
end

% QPSK
M = 4; k = log2(M);
const_QPSK = exp(1j * (pi/4 + 2*pi*(0:M-1)/M));  

for i = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(i)/10);
    EsN0 = EbN0 * k;
    noise_var = 1 / (2 * EsN0);

    err = 0; iter = 0;
    while err < target_err && iter < iter_max
        bits = randi([0 1], 1, num_bits);
        pad = mod(length(bits), k);
        if pad ~= 0
            bits = [bits, zeros(1, k - pad)];
        end
        bits_tx = bits;
        symbols = reshape(bits, k, []).';
        idx = bi2de(symbols, 'left-msb') + 1;
        s = const_QPSK(idx);
        noise = sqrt(noise_var) * (randn(1, length(s)) + 1j * randn(1, length(s)));
        r = s + noise;

        d = abs(r.' - const_QPSK);
        [~, idx_hat] = min(d, [], 2);
        bits_hat = de2bi(idx_hat - 1, k, 'left-msb');
        bits_hat = reshape(bits_hat.', 1, []);
        err = err + sum(bits_tx ~= bits_hat(1:length(bits_tx)));
        iter = iter + 1;
    end
    BER_QPSK(i) = err / (length(bits_tx) * iter);
end

% 8PSK
M = 8; k = log2(M);
const_8PSK = exp(1j * 2 * pi * (0:M-1)/M);

for i = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(i)/10);
    EsN0 = EbN0 * k;
    noise_var = 1 / (2 * EsN0);

    err = 0; iter = 0;
    while err < target_err && iter < iter_max
        bits = randi([0 1], 1, num_bits);
        pad = mod(length(bits), k);
        if pad ~= 0
            bits = [bits, zeros(1, k - pad)];
        end
        bits_tx = bits;
        symbols = reshape(bits, k, []).';
        idx = bi2de(symbols, 'left-msb') + 1;
        s = const_8PSK(idx);
        noise = sqrt(noise_var) * (randn(1, length(s)) + 1j * randn(1, length(s)));
        r = s + noise;

        d = abs(r.' - const_8PSK);
        [~, idx_hat] = min(d, [], 2);
        bits_hat = de2bi(idx_hat - 1, k, 'left-msb');
        bits_hat = reshape(bits_hat.', 1, []);
        err = err + sum(bits_tx ~= bits_hat(1:length(bits_tx)));
        iter = iter + 1;
    end
    BER_8PSK(i) = err / (length(bits_tx) * iter);
end


%plots

semilogy(EbN0_dB, BER_BPSK, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BPSK'); hold on;
semilogy(EbN0_dB, BER_QPSK, 'r-', 'LineWidth', 1.5, 'DisplayName', 'QPSK');
semilogy(EbN0_dB, BER_8PSK, 'g-', 'LineWidth', 1.5, 'DisplayName', '8PSK');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend('Location','southwest');
title('BER BPSK, QPSK y 8PSK');
