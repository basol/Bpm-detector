function BPMDetection()
tic;
clc,clear,close all;

[file, path] = uigetfile('F:\Matlab works\*.wav', 'Select a .wav file');
if isequal(file,0)
   disp('User selected Cancel');
   return;
else
   disp(['User selected ', fullfile(path,file)]);
end


SampleAudio = fullfile(path,file);
BPM = Main(SampleAudio)
toc;
end
	
%% Main
function [BPM] = Main(SampleAudio)

[AudioSignal, SampleRate] = audioread(SampleAudio);

% Selection first channel of audio stereo -> mono
AudioSignal = AudioSignal(:,1);


%% Initialize parameters
MinBPM = 60;                % Minimum BPM
MaxBPM = 250;               % Maximum BPM
EnvelopeReducer = 250;      % Envelope reducer
NumSubBand = 6;             % Number of Subband

%% Scheirer octave ranges
SRTP = [1,200,400,800,1600,3200];
ENDP = [200,400,800,1600,3200,6400]; 
Spacing = round((SampleRate/2)/EnvelopeReducer);

%%
SumEnvelope = 0;
for i = 1:NumSubBand
	
	% Calculation subband
	[SubBand] = SubBandFunction(AudioSignal,SampleRate,SRTP(i),ENDP(i));

	% Envelop construction
	Envelope = EnvelopeFunction(SubBand, i, SampleRate);

	% SampleRate Reduction
	EnvelopeReduced = Envelope(1:Spacing:length(Envelope)); % > 100Hz 

	% Sum envelope
	SumEnvelope = SumEnvelope + EnvelopeReduced;

end

%% Autocorrelation on SumEnvelope
ACRL = AutoCorrelation(SumEnvelope, EnvelopeReducer, MinBPM, MaxBPM);

% Peak detection
[~, Index] = max(ACRL);

% Calculation BPM
BPM = 60*EnvelopeReducer/(Index);

end

%% Calculation subband
function [SubBand] = SubBandFunction(AudioSignal, SampleRate, L, H)

% Low frequency band / Cutoff frequency / cutoffFreq = Freq / (samplerate/2)
LFBand = L/(SampleRate/2); 

% Creating Butterworth filter
[b, a] = butter(2, LFBand, 'high');

% Applying filter
FilteredAudioSignal = filtfilt(b, a, AudioSignal);

% High Frequency Band / Cutoff frequency
HFBand = H/(SampleRate/2);

% Creating Butterworth filter
[b, a] = butter(2, HFBand, 'low');

% Applying filter
SubBand = filtfilt(b, a, FilteredAudioSignal);

end

% Envelop Construction
function [Envelope] = EnvelopeFunction(SubBand, i, SampleRate)

% Full Wave Retification / rectify-and-smooth method 
FWR = abs(SubBand);

% Downsampling
DownSampling = downsample(FWR, 2);

% Normalization (mean removal)
Normalization = DownSampling - mean(DownSampling);

% Eliminate the noise in the subband
Nw=0.1*SampleRate/2;
w=ones(Nw,1)/Nw;     
Envelope = conv(Normalization, w,'same');

end

%% Autocorrelation on Envelope
function [ACRL] = AutoCorrelation(Envelope, EnvelopeReducer, MinBPM, MaxBPM)

StartP = round((60 * EnvelopeReducer) / (MaxBPM));
EndP = round((60 * EnvelopeReducer) / (MinBPM));
SampleNum = length(Envelope) - EndP;

ACRL = zeros(EndP,1);
for k = StartP:EndP
	Sum = 0;
	for i = 1:SampleNum
		Sum = Sum + (Envelope(i)*Envelope(i + k));
	end
	ACRL(k) = Sum;
end

end

