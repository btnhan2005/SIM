function [sig_out] = MMSE_Equalization(sig_in, H, SNR_dB, Pt)
%MMSE_EQUALIZATION
SNR = 10^(SNR_db/10);
N_0 = Pt/SNR;
I = eye(sig_in);
G = H' * inv(H' * H + N_0 * I);
sig_out = G * sig_in;
end

