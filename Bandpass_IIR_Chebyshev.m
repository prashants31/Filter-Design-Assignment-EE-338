%--------------------------------------------------------------------------------------------------------------------------------
%Prashant Shettigar 18D070063
%EE338 Digital Signal Processing
%Filter Design Assignment
%IIR Bandpass Filter
%--------------------------------------------------------------------------------------------------------------------------------

%Band Edge speifications
F_sl = 61.7e3;
F_pl = 65.7e3;
F_ph = 85.7e3;
F_sh = 89.7e3;

%Transformed Band Edge specs using Bilinear Transformation
F_samp = 330e3;                 %Sampling frequency
Ws1 = tan(F_sl/F_samp*pi);          
Wp1 = tan(F_pl/F_samp*pi);
Wp2 = tan(F_ph/F_samp*pi);
Ws2 = tan(F_sh/F_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(Wp1*Wp2);
B = Wp2-Wp1;

%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+1i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-1i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+1i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-1i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5



%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coefficients of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coefficients of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,2048*2048, F_samp);
plot(f,abs(H))

%plotting horizontal and vertical lines corresponding to certain frequency
%and magnitude response
yline(0.15,'-r','Magnitude = 0.15')
yline(0.85,'-r','Magnitude = 0.85')
yline(1,'-r','Magnitude = 1')
xline(61.7e3,'-m','f = 61.7kHz')
xline(65.7e3,'-m','f = 65.7kHz')
xline(85.7e3,'-m','f = 85.7kHz')
xline(89.7e3,'-m','f = 89.7kHz')

xlim([20000,130000])
ylim([0,1.5])
xlabel('Frequency (in Hz)')
ylabel('Magnitude')
title('Magnitude Plot')
hold on
grid on