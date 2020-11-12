%--------------------------------------------------------------------------------------------------------------------------------
%Prashant Shettigar 18D070063
%EE338 Digital Signal Processing
%Filter Design Assignment
%FIR Bandpass Filter
%--------------------------------------------------------------------------------------------------------------------------------

%Sampling frequency
F_samp = 330e3;

%Band Edge speifications
F_sl = 61.7e3;
F_pl = 65.7e3;
F_ph = 85.7e3;
F_sh = 89.7e3;

%normalized edge frequency
w_pl = F_pl*2*pi/F_samp;
w_ph = F_ph*2*pi/F_samp;
w_sl = F_sl*2*pi/F_samp;
w_sh = F_sh*2*pi/F_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*0.0242*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 45;

%Ideal bandpass impulse response of length "n"
w_c1 = (w_pl + w_sl)/2;
w_c2 = (w_ph + w_sh)/2;
bp_ideal = ideal_lp(w_c2,n) - ideal_lp(w_c1,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, F_samp);
plot(f,abs(H))

%plotting horizontal and vertical lines corresponding to certain frequency
%and magnitude response
yline(0.15,'-r','Magnitude = 0.15');
yline(0.85,'-r','Magnitude = 0.85');
yline(1.15,'-r','Magnitude = 1.15');
xline(61.7e3,'-m','f = 61.7kHz');
xline(65.7e3,'-m','f = 65.7kHz');
xline(85.7e3,'-m','f = 85.7kHz');
xline(89.7e3,'-m','f = 89.7kHz');

xlim([20000,130000])
ylim([0,1.5])
xlabel('Frequency (in Hz)')
ylabel('Magnitude')
title('Magnitude Plot')
hold on
grid on 