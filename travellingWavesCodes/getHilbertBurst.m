% Modified from getBurstLengthHilbert in https://github.com/supratimray/GammaLengthProjectCodes.

function [burstLengthS,burstStartS,burstEndS,burstTS,bpfSignal,hilbertPower, thresholdFactor] = getHilbertBurst(analogData,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,burstFreqRangeHz,filterOrder,req)

if ~exist('stimulusPeriodS','var');     stimulusPeriodS=[0.5 1.5];      end
if ~exist('baselinePeriodS','var');     baselinePeriodS=[-1 0];         end
if ~exist('burstFreqRangeHz','var');    burstFreqRangeHz=[40 60];       end

numTrials=size(analogData,1);

hilbertPower=zeros(size(analogData));
bpfSignal=zeros(size(analogData));
for i=1:numTrials
    [hilbertPower(i,:),bpfSignal(i,:)] = getHilbertPower(analogData(i,:),timeVals,burstFreqRangeHz,filterOrder,req);
end

mBL = getMedianBaseline(hilbertPower,timeVals,baselinePeriodS);
if isempty(thresholdFactor)
    mST = getMedianBaseline(hilbertPower,timeVals,stimulusPeriodS);
    thresholdFactor = mean(mST./mBL);
    disp(['Using threshold factor of: ' num2str(thresholdFactor)]);
end

threshold=thresholdFactor*mBL;

burstLengthS=cell(1,numTrials);
burstStartS=cell(1,numTrials);
burstEndS=cell(1,numTrials);

for i=1:numTrials
    if displayFlag
        clf;
        disp(['Trial : ' num2str(i) ' of ' num2str(numTrials)]);
        subplot(211);plot(timeVals,bpfSignal(i,:));xlim([stimulusPeriodS(1)-0.5 stimulusPeriodS(2)+0.5]);
    end
    [burstLengthS{i},burstStartS{i},burstEndS{i},burstTS(i,:)] = getBurstLengthHilbertSingleTrial(hilbertPower(i,:),timeVals,threshold,stimulusPeriodS,displayFlag);
end
end

function [burstLengthS,burstStartS,burstEndS,burstTS] = getBurstLengthHilbertSingleTrial(bandPassPowerSingleTrial,timeVals,threshold,stimulusPeriodS,displayFlag)
stPos = intersect(find(timeVals>=stimulusPeriodS(1)),find(timeVals<stimulusPeriodS(2)));
st=timeVals(1,stPos);
stBandPassPowerSingleTrial=bandPassPowerSingleTrial(1,stPos);

seedPosList=find((stBandPassPowerSingleTrial/threshold)>1); % Seed point is any time point which exceeds the threshold
burstStartPosList = [];
burstEndPosList = [];
timePos=1;
timeLength=length(st);
while (timePos < timeLength) % As long as the time position is less than the search interval
    tPos = seedPosList(find(seedPosList>timePos,1)); % Find the position of the next seed
    minPos = find((stBandPassPowerSingleTrial(1:tPos)/(threshold/2))<1,1,'last');
    if ~isempty(tPos)
        maxPos = find((stBandPassPowerSingleTrial(tPos+1:timeLength)/(threshold/2))<1,1,'first');

        if isempty(minPos)
            tMinPos=1;
        else
            tMinPos=minPos;
        end
        if isempty(maxPos)
            tMaxPos=timeLength;
        else
            tMaxPos=tPos+maxPos;
        end
        timePos=tMaxPos;
        burstStartPosList=cat(1,burstStartPosList,stPos(tMinPos));
        burstEndPosList=cat(1,burstEndPosList,stPos(tMaxPos));
    else
        timePos=timeLength;
    end
end
% 
burstTS = nan(1,length(bandPassPowerSingleTrial));
for ind = 1:length(burstStartPosList)
    burstTS(burstStartPosList(ind):burstEndPosList(ind)) = 1;
end
burstStartS = timeVals(burstStartPosList);
burstEndS = timeVals(burstEndPosList);
burstLengthS = burstEndS-burstStartS;





if displayFlag==1
    subplot(212);plot(timeVals,log10(bandPassPowerSingleTrial));hold on;
    line([timeVals(1),timeVals(length(timeVals))],[log10(threshold),log10(threshold)],'color','r');
    line([timeVals(1),timeVals(length(timeVals))],[log10(threshold/2),log10(threshold/2)],'color','g');
    for j=1:length(burstStartPosList)
        line([timeVals(burstStartPosList(j)),timeVals(burstEndPosList(j))],[log10(threshold/2),log10(threshold/2)],'lineWidth',3,'color','k');
        plot(timeVals(burstStartPosList(j)),log10(threshold/2),'square','linewidth',1,'MarkerSize',8,'color','k');
    end
    xlim([stimulusPeriodS(1)-0.5 stimulusPeriodS(2)+0.5]);
    pause;
    clf;
end
end

function mBL = getMedianBaseline(smoothedPowerBpfSignal,timeVals,baselinePeriodS)
blPos = intersect(find(timeVals>=baselinePeriodS(1)),find(timeVals<baselinePeriodS(2)));
blData = smoothedPowerBpfSignal(:,blPos);
mBL = median(blData(:));
end

function [hilbertPower,bpfSignal] = getHilbertPower(signal,timeVals,gammaFreqRangeHz,filtOrder,req)
% This function computes the instantaneous power at gamma band using
% Hilbert transform. It uses a Butterworth filter to bandpass filter the
% signal in the gamma range

% Input - signal              : Signal in one row
%         timeVals            : Time Values in a row (seconds)
%         gammaFreqRangeHz    : Frequency range for which the Wavelet
%                               transform should be computed
%         filtOrder           : Order of the Butterworth filter
%         Req                 : If Req is 1, data is filtered, else no filtering is done 
% Output - hilbertPower       : Instantaneous power calculated using
%                               Hilbert transform
%          bpfSignal          : Signal bandpass filtered in gamma band
%                               using Butterworth filter

%BPF the signal
if req==1
    Fs=1/(timeVals(2)-timeVals(1));
    normBand=gammaFreqRangeHz/(Fs/2);
    [b,a]=butter(filtOrder,normBand,'bandpass');
    bpfSignal=filtfilt(b,a,signal);
else 
    bpfSignal = signal;
end
%Find Inst. power using HT
hilbSignal=hilbert(bpfSignal);
hilbertPower=abs(hilbSignal).^2;
end