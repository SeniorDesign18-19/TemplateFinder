clear
clc
load hfo_20_awake.mat

% Filters EEG data with band-pass filter from 1 Hz to 35 Hz
for i = 1:18
    XNew(:,i) = downsample(X(:,i),4);
    [b,a] = butter(3,[1/256 35/256],'bandpass');
    XNew(:,i) = filter(b,a,XNew(:,i));
    Sig2Chan(i) = var(XNew(:,i));
end

DataStart = 151;
DataEnd = 301;
TotalSampleNumber = length(XNew)/150*3;

index = 1;

tic;
VectorFinal = {};
while DataEnd < length(XNew)
    
    for Channel = 1:18
        
        Energy(Channel,index) = 0;
        LineLength(Channel,index) = 0;
        
        for i = 1:150
            Energy(Channel,index) = XNew(DataStart + i - 1,Channel)^2 + Energy(Channel,index);
            LineLength(Channel,index) = abs(XNew(DataStart + i - 2,Channel) - XNew(DataStart + i - 1,Channel)) + LineLength(Channel,index);
        end
        RMS(Channel,index) = rms(XNew(DataStart:DataEnd,Channel));
        ZeroCrossingRate(Channel,index) = ZCR(XNew(DataStart:DataEnd,Channel));
        
        SamplePrev = XNew((DataStart-50):(DataEnd-50),Channel); % Previous window
        Sig2Sam(Channel,index) = var(X(DataStart:DataEnd,Channel)); % Variance of sample
        Sig2SamPrev(Channel,index) = var(SamplePrev); % Variance of previous sample
        FVarPrev(Channel,index) = Sig2Sam(Channel)/(Sig2Sam(Channel) + Sig2SamPrev(Channel,index)); % Variance ratio between sample and previous sample
        FVarChan(Channel,index) = Sig2Sam(Channel)/(Sig2Sam(Channel) + Sig2Chan(Channel)); % Variance ratio between channel and sample

        AbsolutePeak(Channel,index) = max(XNew(DataStart:DataEnd,Channel));
        PeakToPeak(Channel, index) = max(XNew(DataStart:DataEnd,Channel)) - min(XNew(DataStart:DataEnd,Channel));
        peakthreshold = 0.1*AbsolutePeak(Channel,index);
        
        pxx = pwelch(XNew(DataStart:DataEnd));
        MeanPxx = mean(pxx);
        PSD(Channel,index) = MeanPxx;
        
        Vector1 = [RMS(Channel,index),ZeroCrossingRate(Channel,index),Sig2Sam(Channel,index),Sig2SamPrev(Channel,index),FVarPrev(Channel,index),FVarChan(Channel,index),AbsolutePeak(Channel,index),PeakToPeak(Channel,index),PSD(Channel,index)];
        VectorFinal{Channel,index} = Vector1;
%         display(Channel);
%         display(index);
    end
    
    DataStart = DataStart + 50;
    DataEnd = DataEnd + 50;
    index = index + 1;
    
end

toc;