clear all
close all
clc

L = 500; 
r = 500; 
n = ceil(0.98 * L); 
m = L-n;

t = linspace(0,r,L);

x = (0.8*cos(0.5*pi*t) + 5*cos(0.2*pi*t));
figure; subplot(2,1,1); plot(t,x); title("Input Signal without noise");

p = 2; 
p_ = 2*p;

% add noise
noise_var = 1;
for i = 1:L
    x_n(i) = x(i) + sqrt(noise_var)*randn;
end

subplot(2,1,2); plot(t,x_n); title("Input Signal corrupted with noise");

% find r_x autocorrelation matrix using method given in paper

X = zeros(n,m); % sample correlation matrix

for i=1:n
    for j=1:m
        X(i,j) = x_n(i-1+j-1 + 1);
    end
end

% X = hankel(x_n(1:m),x(m+1:end));

r_x = X' * X;

[V,D] = eig(r_x);

% V_exp = V;
% D_exp = D;

% sort eigen vectors and eigen values in descending order

for i = 1:floor(size(D,1)/2)
    temp = D(i,i);
    D(i,i) = D(size(D,1)-i+1,size(D,1)-i+1);
    D(size(D,1)-i+1,size(D,1)-i+1) = temp;
end

for i = 1:floor(size(V,2)/2)
    temp = V(:,i);
    V(:,i) = V(:,size(V,2)-i+1);
    V(:,size(V,2)-i+1) = temp;
end

omega = 0:0.01:pi;

e = zeros(m,size(omega,1)); % a single column of e corresponds to the steering vector of a given frequency

for c=1:size(omega,2)
    for r_=1:m
        e(r_,c) = real(exp(1i*omega(c)*(r_-1)));
    end
end

d_squared = zeros(1,size(omega,2));

for c=1:size(omega,2)
    steering_vector = e(:,c);
    sum = 0;
    for k=p_+1:m
        % V_ = fftshift(fft(V(:,k)));
        sum = sum + ( abs( steering_vector' * V(:,k) ) ) ^ 2;
    end
    d_squared(c) = sum;
end

pseudo_spect = 1 ./ d_squared;

pseudo_spect_db = 10*log(abs(pseudo_spect));

% figure; subplot(1,2,1); plot(omega/pi,abs(d_squared)); title("d^2"); subplot(1,2,2); plot(omega/pi,abs(pseudo_spect)); title("Pseudo Spectrum");

figure; plot(omega/pi,pseudo_spect_db,'r'); title("Pseudo Spetrum in dB"); hold on;

[~,R_X] = corrmtx(x_n,m-1,'mod');
pmusic(R_X,p_);
legend("Experimental","MATLAB in-built");