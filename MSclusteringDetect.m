%%
% this script shows how to use the microsaccade detection method
% published in Otero-Millan et al. Journal of Vision 2014
% Set up variables --------------------------------------------------------
folder = fileparts(mfilename('fullpath'));
if ( isempty( folder) )
    folder = pwd;
end
folder = [folder '\MSclustering_tmp'];
session = 'test';
samplerate = 1000;

% load('tmpData.mat')
timestamps=1:size(tmpData,2);
trialsidx=1:size(tmpData,1);
tmpData(:,:,5)=repmat(timestamps,size(tmpData,1),1);
tmpData(:,:,6)=repmat(trialsidx',1,size(tmpData,2));

% concatinate trials
tmp2=(permute(tmpData,[3 2 1]));
tmp3=tmp2(:,:);
samplestmp=tmp3(:,isnan(sum(tmp3))==0)';
trialvector=samplestmp(:,6);
%% lowpass filter
windowSize = 12;
eyetmpmat=zeros(size(samplestmp(:,1:4)));

for ie=1:4
    data=samplestmp(:,ie);
    data2=filter(ones(1,windowSize)/windowSize,1,[mean(data)*ones(windowSize,1);data]);
    eyetmpmat(:,ie)=data2((windowSize+1):end);
end

% transfer into degree of visual angle 
scale=100;
samples=[samplestmp(:,5),eyetmpmat./scale];
blinks=zeros(length(samples),1);
%
% Loads the recording and prepares it por processing
recording = ClusterDetection.EyeMovRecording.Create(folder, session, samples, blinks, samplerate);

% Runs the saccade detection
[saccades, stats] = recording.FindSaccades();

% Plots a main sequence
enum = ClusterDetection.SaccadeDetector.GetEnum;
figure
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');
%% plot some ms
timelimt=8;
saccades2=saccades(saccades(:,enum.duration)>timelimt,:);
MSnum=size(saccades2,1);
trialtmp=trialvector(saccades2(:,enum.startIndex));
Trialabel=unique(trialtmp);
TRafterons=Trialabel;
TRnum=length(Trialabel);
for icount=1:ceil(TRnum/9)
    numT=9*icount;
    scrsz=get(0,'ScreenSize');% get screen size for output display
    
    h1=figure('NumberTitle','off','Name','Check Trials','Position',scrsz);
    for ip=1:9
        iMS=ip+numT;%randi(MSnum);
        
        subplot(3,3,ip)
        iTrial=Trialabel(iMS);
        trialstrart=find(trialvector==iTrial);
        selvect=trialtmp==iTrial;
        eyeMatmp=samples(trialstrart,2:end);
        plot(eyeMatmp);hold on
        selvect2=find(selvect==1);
        selMS=saccades2(selvect2,:);
        for ims=1:size(selMS,1)
        saccstart=selMS(ims,enum.startIndex)-trialstrart(1);
        saccend=selMS(ims,enum.endIndex)-trialstrart(1);
        plot(saccstart:saccend,eyeMatmp(saccstart:saccend,:),'k','linewidth',2)
        end
        if sum(selMS(:,enum.startIndex)>(trialstrart(1)+300))==0
            TRafterons(iMS)=NaN;
        end
        title(['Trial ' num2str(iTrial)])
    end
    % close all
    pause
end
%% orientation
close all

figure;
subplot(3,1,1);hist(saccades2(:,enum.rightDirection),180)
subplot(3,1,2);hist(saccades2(:,enum.leftDirection),180)
subplot(3,1,3);hist(saccades2(:,enum.direction),180)
