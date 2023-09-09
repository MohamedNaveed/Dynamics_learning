function [output] = fft_signal(signal, model, fig_title)

T = model.dt; %time period
Fs = 1/T; %sampling frequency

L = length(signal)+1;
%NFFT = L;%1024;

fft_signal = fft(signal, L, 2);

%fft_signal_shifted = fftshift(fft_signal);
%fVals = Fs*(-NFFT/2:NFFT/2-1)/NFFT;

P2 = abs(fft_signal/L);
P1 = P2(:,1:L/2 + 1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
f = Fs*(0:(L/2))/L;

temp = phase(fft_signal);
signal_phase = temp(:,1:L/2 + 1);

if size(signal,1) == 2
    figure;
    subplot(2,1,1)
    plot(f,abs(P1(1,:)),'LineWidth',3);
    ylabel('theta');
    title(fig_title)
    subplot(2,1,2)
    plot(f,abs(P1(2,:)),'LineWidth',3);
    ylabel('theta dot');
    xlabel('freq');
else
    
    figure;
    %subplot(2,1,1)
    plot(f,P1(1,:),'LineWidth',3);
    ylabel('magnitude');
    title(fig_title)
    xlabel('freq');
    %subplot(2,1,2)
    %plot(f,signal_phase,'LineWidth',3);
    %ylabel('Phase');
    %xlabel('freq');
    
end



end

