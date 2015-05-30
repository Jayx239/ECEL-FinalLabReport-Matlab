%%
% Jason Gallagher
% ECEL-301 section 066
% Lab 8
% Dr. Gail Rosen
% 5/28/2015
clc,clear all;

%%

%% Calculate Component values
% Specification values set 2
AvVal = -10; % Pass Band Gain Value
R1 = 2400;  % input Resistor
r=8100;     % -3dB frequency
% --------------------------
%Calculated Values
Rf = -AvVal*R1;    % Feedback resistor calculation base on Av = -RF/R1
% Ideal Capacitance Calculation
C = 1/(2*pi*r*Rf);

%% Calculate Ideal Gain Vs Frequency Curve
% Import Excel Data
ExGainVV = xlsread('LabViewData','LabViewData','E3:E102');
ExFreq = xlsread('LabViewData','LabViewData','A3:A102');
% Calculate Ideal Gain Values
freq = logspace(1,6,100);
Av = 20*log10((Rf/R1)./(sqrt(1+(2*pi.*freq*Rf*C).^2)));
semilogx(freq,Av,ExFreq,ExGainVV);  % Plot Gain vs Freq

hold on;
%% Get Cutoff Frequencys
% Theoretical
CutoffFreqTheo = 0; % Store Cutoff frequency
indexTheo = 0;  % index for cutoff used to find gain at frequency
for i=1:100
    if Av(i) <= 17
       CutoffFreqTheo = freq(i);    % get closest frequency to -3db
       indexTheo = i;               % get index
       break;                       % Stop loop
    end;
end;
% Display Values
display('Theoretical Cutoff Frequency');
display(CutoffFreqTheo);
display('Theoretical Cutoff Gain Value');
display(Av(indexTheo));
% Add cutoff point to plot
plot(CutoffFreqTheo,Av(indexTheo),'bo');
% Experimental
CutoffFreqEx = 0;
indexEx = 0;
for i=1:100
    if ExGainVV(i) <= 17
       CutoffFreqEx = ExFreq(i);
       indexEx = i;
       break; 
    end;
end;
% Display Values
display('Experimental Cutoff Frequency');
display(CutoffFreqEx);
display('Experimental Cutoff Gain Value');
display(ExGainVV(indexEx));
% add cutoff point to plot
plot(CutoffFreqEx,ExGainVV(i),'go');

plot(freq(1:indexTheo-1),Av(1:indexTheo-1),'ko');   % Theo Pass Band
plot(ExFreq(1:indexEx-1),ExGainVV(1:indexEx-1),'ro');  % Exp Pass Band
% Add Labels
title('Gain vs Frequency');
format('short')
legend('Theoretical Gain', 'Experimental Gain',...
    strcat('Theoretical Corner Frequency (',num2str(CutoffFreqTheo),' Hz)'),...
    strcat('Experimental Corner Frequency (',num2str(CutoffFreqEx),' Hz)'),...
    'Pass Band Theoretical','Pass Band Experimental',...
    'Location', 'SouthWest');
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');

%% Pass Band Gain
% Gather all points prior to cutoff 
passBandTheoVals = Av(1:indexTheo); % All points before corner freq
passBandExVals = ExGainVV(1:indexEx);   % All points before corner freq
passGainTheo = mean(passBandTheoVals);  % calculated using mean value
passGainEx = mean(passBandExVals);  % Calculated using mean
display(passGainTheo);  % Display
display(passGainEx);    % Display

%% Calculate Slopes of stopband (Roll off)
%Theoretical
figure(2);
stopBandTheoVals = rot90(Av((indexTheo):100));      % Gain Vals
stopBandFreqTheoVals = rot90((freq((indexTheo):100))) %Frequency Vals
LinearFitTheo = fit(stopBandFreqTheoVals,stopBandTheoVals,'poly1');
plot(LinearFitTheo,'b');
hold on;
%
display('Theoretical Stop Band Gain');
display(LinearFitTheo(1));


% Experimental
stopBandExVals = (ExGainVV(indexEx:100));   % Gain Vals
stopBandFreqExVals = (ExFreq(indexEx:100)); % Frequency vals
LinearFitEx = fit(stopBandFreqExVals,stopBandExVals,'poly1');
plot(LinearFitEx,'g');
display('Experimental Stop band Gain');
display(LinearFitEx(1));
plot(freq((indexTheo):100),Av((indexTheo):100),'b*'); % Orig theo roll off
plot(stopBandFreqExVals,stopBandExVals,'g*');   % Orig Ex roll off
title('Roll Off Gain Vs Frequency');
legend(strcat('Theoretical Gain Fit (',num2str(LinearFitTheo(1)),' V/V)'),...
    strcat('Experimental Gain Fit (',num2str(LinearFitEx(1)),' V/V)')...
    ,'Theoretical Gain', 'Experimental Gain');

xlabel('Frequency (Hz)');
ylabel('Gain (dB)');

%%
%%stopBandTheoValsAdj = rot90(Av((indexTheo-12):100));      % Gain Vals adjusted
%%stopBandFreqTheoValsAdj = rot90((freq((indexTheo-12):100))) %Frequency Vals adjusted
%%LinearFitTheoAdj = fit(stopBandFreqTheoValsAdj,stopBandTheoValsAdj,'poly1');
%%plot(LinearFitTheoAdj,'r');
%%hold on;
%

%%display('Adjusted Theoretical Stop Band Gain');
%%display(LinearFitTheoAdj(1));
%

%strcat('Theoretical Gain Adjusted (',num2str(LinearFitTheoAdj(1)),' V/V)')

%% Resources used
% http://en.wikipedia.org/wiki/Low-pass_filter