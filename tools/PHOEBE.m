function [keep,SCI,xpower] = PHOEBE(din, cutoff, thres1, thres2, wlens)
%This calculates two thresholds to assess signal quality of fNIRS data.
%Similar as in Pollonini et al 2016: https://doi.org/10.1364/BOE.7.005104
%But here, I considered a shifting window of 10s for the computation.
%Note: original assessment was suggested for two wavelengths only.
%
%Last Update: 2018-03-04
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

if nargin < 2
    cutoff = [0.5 2.5]; %retrieve cardiac frequency
end

if nargin < 3
    thres1 = 0.8; %cross-correlation at 0 lag
end

if nargin < 4
    thres2 = 0.1; %power spectrum peak
end

if nargin < 5
    wlens = 10; %10s sliding window
end

%retrieve #sessions and #channels
%note: assumed #channels constant
subjs = length(din);
chtot = size(din(1).data,2)/2;

%initialize variables
keep = zeros(subjs,chtot);
SCI = zeros(subjs,chtot);
xpower = zeros(subjs,chtot);

%make copy of input data (for filtering purposes)
d = din;

%iterate over sessions provided
for n=1:subjs
    
    %sampling frequency
    Fs = din(n).Fs;
    
    %convert window length to data points
    wlen = round(wlens*Fs);
    
    %set and apply 4th-order Butterworth filter
    order = 4;
    cutoff(2) = min(cutoff(2),(Fs/2)-0.1); %max cutoff < Fs/2 (Nyquist)
    [a,b] = butter(order/2,cutoff./(Fs/2));
    d(n).data = filtfilt(a,b,din(n).data);
    
    %retrieve indexes for each wavelength
    wl_type = unique(din(n).probe.link.type);
    wl1 = d(n).data(:,din(n).probe.link.type==wl_type(1));
    wl2 = d(n).data(:,din(n).probe.link.type==wl_type(2));
    
    %initialize SCI and xpower
    SCI_w = zeros(1,length(1:wlen:(size(wl1,1)-wlen)));
    xpower_w = SCI_w;
    
    %iterate over all channels
    for ch=1:chtot
        
        %initialize counter for sliding window iteration
        i=1;
        
        %run sliding window with length wlen
        for w=1:wlen:(size(wl1,1)-wlen)
            
            %normalize data of each wavelength to its standard deviation
            wl1_norm = wl1(w:w+wlen,ch)./std(wl1(w:w+wlen,ch),1);
            wl2_norm = wl2(w:w+wlen,ch)./std(wl2(w:w+wlen,ch),1);
            
            %calculate scalp coupling index
            SCI_w(i) = xcorr(wl1_norm,wl2_norm,0,'coeff');
            
            xcorr_wl = xcorr(wl1_norm,wl2_norm,'coeff');
            
            %Compute Power Spectrum of normalized cross correlation
            X = xcorr_wl;
            t = 1/Fs:1/Fs:length(X)/Fs;
            y = fft(X);
            P2 = abs(y/length(t));
            P1 = P2(1:round(length(t)/2));
            
            %calculate xpower
            xpower_w(i) = max(P1);
            
            %add to counter
            i = i+1;
             
        end
        
        %final SCI and xpower values (median from sliding window)
        SCI(n,ch) = median(SCI_w);
        xpower(n,ch) = median(xpower_w);
        
        %signal quality - check if both are above the set threshold
        if SCI(n,ch) > thres1 && xpower(n,ch) > thres2
            keep(n,ch) = 1;
        end
        
    end
    
end

end