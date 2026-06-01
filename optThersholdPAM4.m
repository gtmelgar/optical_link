function threshold = optThersholdPAM4(downSampledSignal, prbsSignalLogical)

sd00_noise_tot = std(downSampledSignal(prbsSignalLogical == -3));
sd01_noise_tot = std(downSampledSignal(prbsSignalLogical == -1));
sd11_noise_tot = std(downSampledSignal(prbsSignalLogical == 3));
sd10_noise_tot = std(downSampledSignal(prbsSignalLogical == 1));

I00 = mean(downSampledSignal(prbsSignalLogical == -3));
I01 = mean(downSampledSignal(prbsSignalLogical == -1));
I11 = mean(downSampledSignal(prbsSignalLogical == 3));
I10 = mean(downSampledSignal(prbsSignalLogical == 1));

threshold(1) = (I00.*sd00_noise_tot+I01.*sd01_noise_tot)./ ...
        (sd00_noise_tot+sd01_noise_tot);
threshold(2) = (I01.*sd01_noise_tot+I11.*sd11_noise_tot)./ ...
        (sd01_noise_tot+sd01_noise_tot);
threshold(3) = (I10.*sd10_noise_tot+I11.*sd11_noise_tot)./ ...
        (sd10_noise_tot+sd11_noise_tot);
