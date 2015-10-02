%% adoptee face old new task
% % screen = 39 cm width, REs. 800*600, 60 Hz
% range EyeLink=32deg, at 70cm visual angle=31.13deg.
% 25pix/deg
% 2pix/mm

% image size 382*390

% learning - 14 faces, recognition -  28 faces
% block order: WC,EA,WC,EA

%% load data
% eyedata [situation echantillon list(i,1) randx randy time eyex eyey diffx diffy ];
% bahvioral [list(i,1) answer timeanswer ctRect(1) ctRect(2) ctRect(3) ctRect(4)];
% situation: WC 1 2 5 6 EA 3 4 7 8| learning 1 3 5 7 recognition 2 4 6 8

subjname='Savina Seo';

%
datafile=[subjname '.txt'];
datatxt=fopen(datafile);
data=textscan(datatxt, '%n%n%n%n%n%n%n%n%n%n');
situation=(data{1});
echantillon=(data{2});
item=(data{3});
randx=(data{4});
randy=(data{5});
timetrial=(data{6});
eyex=(data{7});
eyey=(data{8});
diffx=(data{9});
diffy=(data{10});

eyeraw=[eyex,eyey];
coordina=[randx,randy];
eyet=eyeraw-coordina;
ySize=382;
xSize=390;
% % check data
%
% rawmap = zeros(ySize, xSize);
coordX = (eyet(:, 2));
coordY = (eyet(:, 1));
indx1=coordX>0 & coordY>0 & coordX<ySize & coordY<xSize;
sampDist=[0;sqrt((coordX(2:end)-coordX(1:end-1)).^2+(coordY(2:end)-coordY(1:end-1)).^2)];
veloctmp=smooth(sampDist(indx1));
threltmp=prctile(veloctmp,80);
% indxtf=sub2ind([ySize, xSize],coordX(indx1),coordY(indx1)); % index each fixation location
% unindx = unique(indxtf);% find unique fixation
% [cotind,whe] = histc(indxtf,unindx); % cumulate fixation with same coordinates.
% rawmap(unindx)=cotind;% number of fixation
% figure;imagesc(rawmap)
summary=[];
for ist=1:8
    isitem=unique(item(situation==ist));
    for it=1:length(isitem) %should be 14 for odd situation and 28 for even situation
        %%
        eyetmp=eyet((situation==ist & item==isitem(it)),:);
        timetmp=timetrial(situation==ist & item==isitem(it));
        % invtmp=[0;timetmp(2:end)-timetmp(1:end-1)];
        coordXt = eyetmp(:, 1);
        coordYt = eyetmp(:, 2);
        indx1t=coordXt>0 & coordYt>0 & coordXt<xSize & coordYt<ySize;
        
        sampdist=[0;sqrt((coordXt(2:end)-coordXt(1:end-1)).^2+(coordYt(2:end)-coordYt(1:end-1)).^2)];
        
        veloctmp=smooth(sampdist(indx1t),10);
        coX2=coordXt(indx1t);coY2=coordYt(indx1t);dur2=timetmp(indx1t);
        fixevtmp=bwlabel(veloctmp<threltmp);
        fixevtmp2=zeros(size(fixevtmp));
        for ifx=1:max(fixevtmp)
            fixk2=find(fixevtmp==ifx);
            if length(fixk2)>25 % long fixations are good to go
                fixevtmp2(fixk2)=fixevtmp2(fixk2)+1;
            elseif ifx==1
                fixk3=find(fixevtmp==(ifx+1));
                if (fixk3(1)-fixk2(end))<5
                    fixevtmp2(fixk2(1):fixk3(end))=fixevtmp2(fixk2(1):fixk3(end))+1;
                end
            elseif ifx==max(fixevtmp)
                fixk1=find(fixevtmp==(ifx-1));
                if (fixk2(1)-fixk1(end))<5
                    fixevtmp2(fixk1(1):fixk2(end))=fixevtmp2(fixk1(1):fixk2(end))+1;
                end
            else
                fixk1=find(fixevtmp==(ifx-1));
                if (fixk2(1)-fixk1(end))<5
                    fixevtmp2(fixk1(1):fixk2(end))=fixevtmp2(fixk1(1):fixk2(end))+1;
                end
                fixk3=find(fixevtmp==(ifx+1));
                if (fixk3(1)-fixk2(end))<5
                    fixevtmp2(fixk2(1):fixk3(end))=fixevtmp2(fixk2(1):fixk3(end))+1;
                end
            end
        end
        fixev=bwlabel(logical(fixevtmp2));
        for ifx=1:max(fixev)
            if sum(fixev==ifx)<25 % short fixations are gone
                fixev(fixev==ifx)=0;
            end
        end
        fixev=bwlabel(logical(fixev));
        
        figure;plot(dur2,veloctmp*10,'k');hold on;
        plot(dur2,coX2,'r');plot(dur2,coY2,'g');ylim([0 xSize])
        plot(dur2(logical(fixev)),coX2(logical(fixev)),'o','markers',2);
        plot(dur2(logical(fixev)),coY2(logical(fixev)),'o','markers',2);
        
        summarytmp=zeros(max(fixev),7);
        for ifx=1:max(fixev)
            fixk=find(fixev==ifx);
            gm=gmdistribution.fit([coX2(fixk) coY2(fixk)],1);
            summarytmp(ifx,[1 2])=gm.mu;
            summarytmp(ifx,3)=dur2(fixk(end))-dur2(fixk(1));% duration
            summarytmp(ifx,4)=dur2(fixk(1));% start
            summarytmp(ifx,5)=ifx;% fixation number
            summarytmp(ifx,6)=it;% item number
            summarytmp(ifx,7)=ist;% block number
        end
        summary=[summary;summarytmp];
    end
end
