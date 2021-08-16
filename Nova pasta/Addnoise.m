function sig_out = Addnoise(sig,order,noise_percent)

sig_out = sig+noise_percent/100*order*randn(length(sig),1); % or sig_out = awgn(sig,40,'measured');
