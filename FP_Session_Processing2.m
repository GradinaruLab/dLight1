function y = FP_Session_Processing2(x,t1,t2)
% Process Fiber Photometry (FP) data according to Lerner et al. Cell, 2015.
% dF/F = (490nm signal - fitted 400nm signal)/fitted 400nm signal
%
% Inputs
% x: FP data file loaded by "TDT2mat.m" (required)
% t1: start of data section in seconds (optional)
% t2: end of data section in seconds (optional)
% If there's no t1 and t2, t1 = 3, t2 = end
% 
% Outputs
% Data = Raw FP signal (490/400)
% data = Lowpass filtered FP signal at "LFcut" Hz (see "Basic Variables")
% fit400 = Fitted 400nm signal using linear transformation
% t = Time vector (t1 = 0, t2 = end)
% dF = Processed, dF/F signal (in %)
% TTL = TTL input (onset) if there is any 
% (I have configured digital input as "DIn4")
%
% Version 2 - Nov 7th, 2017
% Original algorithm in Lerner et al., Cell 2015
% Coded by Ryan Cho/Alon Greenbaum, Gradinaru Lab, Caltech

narginchk(1,3) % error if # of input is less than 1 and more than 3

%% Basic Variables and Setup
start = 3; % taking account of manually changing 490 channel frequency
% (seriously, why does it start at 2 Hz at the beginning of recording?)
fs = 382; % sampling rate (in integer)
LFcut = 25; % cut-off frequency of lowpass filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for DA sensor
order = 4; % N-th order for butterworth filter

Data = double(x.streams.LMag.data);
if nargin == 1
    t1 = start; 
    Data = Data(:,t1*fs:end-fs*0.5); 
elseif nargin == 2
    Data = Data(:,t1*fs:end-fs*0.5);
elseif nargin == 3
    Data = Data(:,t1*fs:t2*fs);
end


%% Process FP Data
data = zeros(size(Data));
for i=1:2
    data(i,:) = lowpass(Data(i,:),LFcut,fs,order); % lowpass filter, see below
end
t = 0:1/fs:length(Data)/fs-1/fs;

%% Process TTL input data, if any (caution: I use DIn4 for digital input)
check = isfield(x.epocs,'DIn4'); % check if there's any TTL input
if check == 1
    TTL = x.epocs.DIn4.onset; % extract TTL onset only
    TTL = TTL - t1; % correct for snippet
else
    TTL = []; % otherwise, save it as empty vector
end
% TTL vector is in seconds

%% Fit and Align 405 signal to 490 one <--------------------------------This part is most useful: about least-squares linear fit
% Modeled as 490signal = a * 400signal + b;
% A*theta = B
% where B = [490signal] and
% A = [400signal one-vector];
% (basically, this is like solving linear equation)
B = data(1,:)';
A = [data(2,:)' ones(length(data),1)];
theta = A\B;
fit400 = theta(1)*data(2,:)+theta(2);

% % Or Alternatively, you can use polyfit too
% P = polyfit(data(2,:),data(1,:),1);
% a = P(1); b = P(2);
% fit400 = a.*data(2,:) + b;

dF = 100*detrend((data(1,:)-fit400)./fit400); % Get dF/F in %

% %% Visualize
% % demeaned 490 and 400 signals
% figure; subplot(2,1,1);
% plot(t,data(1,:)-mean(data(1,:))); hold on; plot(t,data(2,:)-mean(data(2,:)),'k');
% set(gca,'FontSize',10); xlim([0 t(end)]);
% legend('GCaMP6 Signal, 490 nm','GCaMP6 Signal, 400 nm');
% xlabel('Time (sec)','FontSize',12);
% ylabel('Fluorescence (au)','FontSize',12);
% title('Lowpass-filtered and de"mean"ed 490 nm and 400 nm signals','FontSize',12);
% % 490 and fit signal
% subplot(2,1,2);
% plot(t,data(1,:)); hold on; plot(t,fit400,'k');
% set(gca,'FontSize',10); xlim([0 t(end)]);
% legend('490 nm','Fitted 400 nm');
% xlabel('Time (sec)','FontSize',12);
% ylabel('Fluorescence (au)','FontSize',12);
% title('490 nm and "fitted" 400 nm signals','FontSize',12);
% 
% % dF/F calcualted using fitted 405 signal
% figure;
% plot(t,dF); hold on;
% set(gca,'FontSize',10); xlim([0 t(end)]);
% xlabel('Time (sec)','FontSize',12);
% ylabel('dF/F (%)','FontSize',12);
% title('Processed GCaMP6 dF/F Signal','FontSize',12);
% % Plot TTL input if it exists
% if isempty(TTL) == 0
%     for i=1:length(TTL);
%         plot([TTL(i) TTL(i)],[min(dF) max(dF)],'r');
%     end
% end
% % If 490 signal is weak - to look like a "flat" line
% if max(dF) <= 3 
%     ylim([-2 5]);
% else
%     ylim([-3 max(dF)+1]);
% end

%% Output
y.Data = Data; % raw data
y.data = data; % lowpassed data
y.fit400 = fit400; % fitted 400nm signal
y.t = t; % time vector
y.dF = dF; % processed dF/F signal
y.TTL = TTL; % external TTL input onset


function y = lowpass(x,cutoff,Fs,n)
% Butterworth Lowpass Filter (zero-phase distortion filter)
% This creates n-th order Butterworth lowpass filter and takes input
% signal and creates output, the filtered signal. 
%
% <Usage>
%
% y = lowpass(x,cutoff,Fs,n)
% 
% where x: the unfiltered, raw signal
%       cutoff: cut-off frequency
%       Fs: sampling rate
%       n: the order of Butterworth filter
%       y: the filtered signal
%
% <Example>
%
% y = lowpass(x,100,2000,4);
%
% Coded by Ryan Cho, Oct 21 2013

[b,a] = butter(n,cutoff/(Fs/2),'low');
y = filtfilt(b,a,x);
return
