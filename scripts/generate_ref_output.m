%% Generate test signals
N = 1000;
k = [0 : 1000];

f1 = 10;
f2 = 20;
f3 = 317;

y1 = sin(2*pi*f1*k/N);
y2 = sin(2*pi*f2*k/N);
y3 = sin(2*pi*f3*k/N);

yt = y1+y2+y3;

%% FFT test output
Y = fft(yt);
%%
fid = fopen('tmp.txt','w');
fprintf(fid, 'const std::vector<std::complex<double>> Yref = {')

for i = 1 : length(Y)
    fprintf(fid, '{%f, %f}', real(Y(i)), imag(Y(i)));
    if i == length(Y)
        fprintf(fid, '};');
    else
        fprintf(fid, ', ');
    end
end
fclose(fid);
%% STFT test signals
% plot(k, yt)

s1 = stft(yt, 'Window', hann(128, 'periodic'), 'OverlapLength', 64, 'FFTLength', 128, 'FrequencyRange', 'onesided').';
s2 = stft(yt, 'Window', hamming(256, 'periodic'), 'OverlapLength', 17, 'FFTLength', 316, 'FrequencyRange', 'onesided').';

%% Spectrogram test signals

b = [1, -0.95];
a = 1;
yt = filter(b, a, yt);

s = spectrogram(yt, hann(128, 'periodic'), 64, 128, 'onesided').';

s = abs(s);

%s = abs(s).^2;
%s = 10 * log(s);

fid = fopen('tmp.txt','w');
fprintf(fid, 'const std::vector<std::vector<double>> yRef = {')
for i = 1 : size(s, 1)
    fprintf(fid, '{');
    for j = 1 : size(s, 2)
        if j == size(s, 2)
            fprintf(fid, '%f', s(i, j));
        else
            fprintf(fid, '%f, ', s(i, j));
        end
    end
    if i == size(s, 1)
        fprintf(fid, '}};');
    else
        fprintf(fid, '},\r\n');
    end
end
fclose(fid);