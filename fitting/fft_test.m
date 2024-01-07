clc;
clear
% Signal parameters
Fs = 1000;          % Sampling frequency
T = 1/Fs;           % Sampling period
L = 1000;           % Signal length
t = (0:L-1)*T;       % Time vector
f = 50;             % Frequency of the signal

% Original signal
x = sin(2*pi*f*t);

% Compute the FFT for the original signal
fftResult = fft(x,4048);

% Compute the power spectrum
powerSpectrum = abs(fftResult/L).^2;

% Plot the power spectrum and mark the peaks
frequencies = Fs*(0:(L/2))/L;
plot(frequencies, powerSpectrum(1:L/2+1));


title('Power Spectrum with Peaks');
xlabel('Frequency (Hz)');
ylabel('Power');
legend('Power Spectrum', 'Peaks');
