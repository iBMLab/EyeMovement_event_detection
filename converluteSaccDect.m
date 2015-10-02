%% set up funtion
Tail=20; % per-saccade/post-saccade intervel, this will the only parameter.
         % it decide the minimal temporal seperation between 2 saccades
xt=-Tail:1:Tail;

% Unit Step Function
yUSF=heaviside(xt)*2-1; % put minial value to -1, now it has range [-1 1]

% High kurtesis function
LPF=@(b,x)(1/b)*exp(-abs(x)/b); % Laplace function
yLPF=LPF(1,xt);
yND= normpdf(xt,0,.45); % Normal distribution (0.45)

% Sigmoid function as saccade model
Sigmoid=@(p,x)p(3)./(1+exp((x+p(1))./p(2)))+p(4);
swift=0; % swift, for deciding the start point of the saccade
ratio=-1; % ratio, for deciding the slop/velocity of the saccade
amplitude=1; % amplitude of the saccade
startpoint=0; % start position of the saccade
p=[swift ratio amplitude startpoint];
ySgM=Sigmoid(p,xt);

% check functions
% figure('NumberTitle','off','Name','check functions','Position',[1 1 1000 200]);
figure;
hold on
plot(xt,yUSF,'k','linewidth',2.5)
plot(xt,yLPF,'b','linewidth',2.5)
plot(xt,yND,'r','linewidth',2.5)
plot(xt,ySgM,'g','linewidth',2.5)
legend({'Unit Step Function';'Laplace Function';'Normal Distribution';'Sigmoid Function'})
%%
eyedata=samplestmp(:,[1 2]);
% eyedata=squeeze(EyeData(1,1:2,:))';

convLayer1=zeros(size(eyedata));
outputLayer1tmp=zeros(size(eyedata));
for ichanal=1:size(eyedata,2)
    convLayer1(:,ichanal)=conv(eyedata(:,ichanal),yUSF,'same');
end
for ieye=1:(size(eyedata,2)/2)
    outputLayer1tmp(:,1+(ieye-1)*2)=sum(abs(convLayer1(:,[1 2]+(ieye-1)*2)),2);
    outputLayer1tmp(:,2+(ieye-1)*2)=sqrt(convLayer1(:,1+(ieye-1)*2).^2+convLayer1(:,2+(ieye-1)*2).^2);
end
if (size(eyedata,2)/2)==1
    outputLayer1=outputLayer1tmp;
else
    outputLayer1(:,1)=mean(outputLayer1tmp(:,[1 3]));
    outputLayer1(:,2)=mean(outputLayer1tmp(:,[2 4]));
end
outputLayer1=outputLayer1(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Peak Prominence method
[pks,loc]=findpeaks(outputLayer1,'MinPeakProminence',200);
% compute the peak prominence and width then kmean clustering with 2
% cluster to identify the saccades

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal envelop method
% y=outputLayer1-repmat(mean(outputLayer1),length(outputLayer1),1);
% outputLayer11=abs(hilbert(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Kmean method
% [kc1,kcm1]=kmeans(outputLayer1,3);
% selchanl=1;
% selkc=find(kcm1(:,selchanl)==min(kcm1(:,selchanl)));
% selkc2=find(kcm1(:,selchanl)==max(kcm1(:,selchanl)));
% outputLayer11=outputLayer1(:,selchanl);
% outputLayer11(outputLayer11<kcm1(1,selchanl))=0;
% outputLayer11(kc1==selkc)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tfce method
% tfce_score(:,1) = tfce2d(outputLayer1(:,1));
% tfce_score(:,2) = tfce2d(outputLayer1(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second convolution method
% The second convolution is not improving the performance, discarde for now
% func=yLPF+abs(cos(xt./2))./10;
% convLayer2=outputLayer1;
% for icycle=1
%     for ichanal=1:size(convLayer2,2)
%         % notmslize
%         %     convLayer2(:,ichanal)=convLayer2(:,ichanal)./max(convLayer2(:,ichanal));
%         convLayer2(:,ichanal)=conv(convLayer2(:,ichanal),func,'same');
%     end
% end
% 
% selchanl=1;
% [kc,kcm]=kmeans(convLayer2,2);
% selkc=find(kcm(:,1)==min(kcm(:,1)));
% outputLayer2=convLayer2(:,selchanl);
% outputLayer2(kc==selkc)=0;
%% ploting
x1=randi(length(eyedata));x1=624123;
% xrange=[x1 x1+3000];

figure('NumberTitle','off','Name','Result');
subplot(3,1,1)
plot(eyedata)
% xlim(xrange)
subplot(3,1,2);
plot(outputLayer1)
hold on
% plot(outputLayer11,'r','linewidth',2.5)
% xlim(xrange)
subplot(3,1,3);
% plot(convLayer2); hold on
% plot(outputLayer2,'r','linewidth',2.5)
% plot(tfce_score)
% plot(diff(outputLayer1))
% xlim(xrange)
%% Saccade Modeling
[pektmp,loc]=findpeaks(outputLayer11);
beta=zeros(size(eyedata,2)/2,length(loc),2,4);
r2=zeros(size(eyedata,2)/2,length(loc),2,2);
beta0=[1 1 1 1];
for ip=1:length(loc)
    tempeak=loc(ip);
    for ieye=1:size(eyedata,2)/2
        for ichannel=1:2
            tmpy=eyedata(xt+tempeak,ichannel+(ieye-1)*2);
            mdl=NonLinearModel.fit(xt',tmpy,Sigmoid,beta0);
            beta(ieye,ip,ichannel,:)=mdl.Coefficients.Estimate;
            r2(ieye,ip,ichannel,:)=[mdl.Rsquared.Adjusted mdl.Rsquared.Ordinary];
        end
    end
end
%%
xt = linspace(min(combin2),max(combin2),1000);
figure;
hist(combin2,xt);
% Get histogram patches
ph = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph);
for i = 1:N_patches
    % Get patch vertices
    vn = get(ph(i),'Vertices');
    % Adjust y location
    vn(:,2) = vn(:,2) + 1;
    % Reset data
    set(ph(i),'Vertices',vn)
end
set(gca,'yscale','log');% Change scale
%% 
% then use boxplot on [combin] to find outliner (threshold)

% then use findpeaks with 'minpeakheight' set as threshold

% for each potential saccade,select the -xtimepoint:xtimpoint from the
% peak,fit model (sigmoid function), the model fitting will the critria
% whether it's a real saccade or just noise. the betas output from fitting
% will the charteristic of the saccade (amplitude, start/end location,
% speed, etc). 
