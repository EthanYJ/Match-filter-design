err_Matched = [];
err_Filtered = [];

for iterator = 1:1:30    
    %% paremeter
    I = iterator*.1; %% noise intensity
    %seq = [-1,1,-1,1,-1,1,1,1,-1,-1,1-1,1-1,1]';
    %% randomly generated training signal
    seq = -5 + (5+5)*rand(45,1);
    for l = 1:length(seq)
        if seq(l)>=0
            seq(l)=1;
        else
            seq(l)=-1;
        end
    end
    barker13 = [1,1,1,1,1,-1,-1,1,1,-1,1,-1,1]';     
    %% generated training signal



    signal = zeros(length(seq),length(barker13))';
    for k = 1:length(seq)
        signal(:,k) =seq(k)*barker13; 
    end
    signal = reshape(signal,1,[]);
    [sig_match,lag] = xcorr(signal,barker13);
    noise=randn(1,length(signal));
    noisy_signal = signal + I*noise;
    [noisy_sig_match,lag_noise_sig] = xcorr(noisy_signal,barker13);

    %plot original signal
%     figure;
%     plot(signal);
%     ylim([-2,2]);
%     title('signal')

    % ploting sig_match
    %figure;
    %plot(lag,sig_match);
    %title('original signal cross_correlation with barker13');
    %plot noisy_signal match
    %figure;
    %plot(lag_noise_sig,noisy_sig_match);
    %title('noisy signal cross_correlation with barker13');

    %% filter designing
    % original training signal correlation Rx
    [Rx,lag_Rx] = xcorr(signal);
    %figure;
    %plot(lag_Rx,Rx);
    %hold on;
    %noise correlation Rv
    [Rv,lag_Rv]= xcorr(noise);
    %plot(lag_Rv,Rv);
    % noise signal correlation Ry
    [Ry,lag_Ry]= xcorr(noisy_signal);
    %plot(lag_Ry,Ry);
    % Rxy
    [Rxy,lag_Rxy]= xcorr(noisy_signal,signal);
    %plot(lag_Rxy,Rxy);
%     hold off;
%     title('signal,noise, noisy signal, signal to noise signal correlation');
%     legend('Rx','Rv','Ry','Rxy');


    %% design a fir filter with order of 64
    l = 64+1;

    %using autocorrelation data matrix size(N+l-1,l) 
    N = length(noisy_signal);
    n= N+l-1;
    %construct training data matrix
    A = zeros(n,l);
    for m = 1:l
        A(m:N+m-1,m)= noisy_signal;
    end

    %Ry matrix
    Ry_mat = zeros(l,l);

    for r = 1:l
       Ry_mat(1:l,r) = Ry(l-r+1:l+l-r);  
    end

    Ry_mat_inv= inv(Ry_mat);
    % here is the filter
    h=Ry_mat\Rxy(l:l+l-1)' ;
    %%%%%%%%%%%%%%%%%%%%%

    % plot filtered signal
    filt_sig = (A(1:N,:)*h)';
    %figure;
    [filt_sig_match,lag_filt_sig]= xcorr(filt_sig,barker13);
    %plot(lag_filt_sig,filt_sig_match);
    %hold on;
%     plot(lag,sig_match);
%     title('training')
    %% test sequence
    % test sequence
    seq_new = [-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,1,-1]';
    signal_new = zeros(length(seq_new),length(barker13))';
    for k = 1:length(seq_new)
        signal_new(:,k) =seq_new(k)*barker13; 
    end
    % new signal
    signal_new = reshape(signal_new,1,[]);

    %new noisy signal
    noise_new=randn(1,length(signal_new));
    noisy_signal_new = signal_new + I*noise_new;

    %new data matrix
    A_new = zeros(n,l);
    for m = 1:l
        A_new(m:N+m-1,m)= noisy_signal_new;
    end

    %new filtered signal
    filt_sig_new = (A_new(1:N,:)*h)';

    % ploting comparing
%     figure;
%     new signal match
    [sig_match_new,lag_new] = xcorr(signal_new,barker13);
    %plot(lag_new,sig_match_new);
%     hold on;
    % new noisy signal match
    [noisy_sig_new_match,lag_noisy_sig_new]= xcorr(noisy_signal_new,barker13);
    %plot(lag_noisy_sig_new,noisy_sig_new_match);

    % new filtered noisy signal match
    [filt_sig_match_new,lag_filt_sig_new]= xcorr(filt_sig_new,barker13);
    %plot(lag_filt_sig_new,filt_sig_match_new);

%     title('testing');
%     legend('original_signal','noisy_signal','filtered_signal');


    %% determine the resulting received sequences
    sentSignal = ((seq+1)./2)';%signal sent
    testSignal = ((seq_new+1)./2)';

    rec_match = [];
    rec_filt_match = [];

    for a = 1:1:length(noisy_sig_match(N-13:end))-1
       if mod(a,13) ~= 0
           continue
       else
           index = a/13;
           if noisy_sig_match(N-13+a) > 0
               rec_match(index) = 1;
           else
               rec_match(index) = 0;
           end
       end
    end

    for a = 1:1:length(noisy_sig_new_match(N-13:end))-1
       if mod(a,13) ~= 0
           continue
       else
           index = a/13;
           if noisy_sig_new_match(N-13+a) > 0
               rec_filt_match(index) = 1;
           else
               rec_filt_match(index) = 0;
           end
       end
    end

    sentSignal;
    rec_match;
    testSignal;
    rec_filt_match;
    %% determine error
    err_Matched(iterator) = sqrt((sentSignal-rec_match)*(sentSignal-rec_match)');
    err_Filtered(iterator) = sqrt((testSignal-rec_filt_match)*(testSignal-rec_filt_match)');
end
%% plot error
figure(1);
plot(.1:.1:3,err_Matched)
hold on
plot(.1:.1:3,err_Filtered)
hold off
title('Error vs Noise Intensity')
xlabel('Noise Intensity (Std. Deviation)')
ylabel('L2 Norm of Error')
%% plot signals
% figure(2);
% plot(sig_match(N-13:end))
% title('Sent Signal')
% figure(3);
% plot(noisy_sig_match(N-13:end))
% title('Received Signal (Matched Filter Only)')
% figure(4);
% plot(noisy_sig_new_match(N-13:end))
% title('Received Signal (Matched Filter + LMS Filter)')

