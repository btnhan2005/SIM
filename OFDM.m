function [sig_out] = OFDM(sig_in, nfft, cplen)
 % dataIn: ma trận (nfft x Nsym), mỗi cột là một OFDM symbol trong miền tần số
    % nfft: kích thước IFFT
    % cplen: độ dài CP
    % IFFT
    x_time = ifft(sig_in, nfft, 1);   % Kích thước: (nfft x Nsym)
    % Thêm CP
    cp = x_time(end-cplen+1:end, :);  % (cplen x Nsym)
    sig_out = [cp; x_time];              % (nfft+cplen x Nsym)
end

