%%
% this script shows how to use the microsaccade detection method
% published in Otero-Millan et al. Journal of Vision 2014

% Set up variables --------------------------------------------------------
folder = fileparts(mfilename('fullpath'));
if ( isempty( folder) )
    folder = pwd;
end
folder = [folder '\data'];
session = 'test';
samplerate = 500;

samples = [];
%  samples(:,1)     timestamps of the recording
%  samples(:,2)     horizontal position of the left eye   
%  samples(:,3)     vertical position of the left eye  
%  samples(:,4)     horizontal position of the right eye     
%  samples(:,5)     vertical position of the right eye  

blinks = [];
%  blinks           binary vector indicating for each sample if it
%                   belons to a blink (1) or not (0)

% Set up variables --------------------------------------------------------


% Loads the recording and prepares it por processing
recording = ClusterDetection.EyeMovRecording.Create(folder, session, samples, blinks, samplerate);

% Runs the saccade detection
[saccades stats] = recording.FindSaccades();

% Plots a main sequence
enum = ClusterDetection.SaccadeDetector.GetEnum;
figure
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');