clc;
clear all;
close all;
% i.
Fs=44.1e3; %define sampling frequency 44.1 kHz
% n=16; %define sample is 16 bits
% ch=1; %define channel mono(1) channel
% a=audiorecorder(Fs,n,ch);%Create object for recording audio; Fs samplerate,n bits,ch channel.
% disp('start spreaking');%This command is start recording
% recordblocking(a,5);%Record a sound with length 5 seconds.
% disp('End of recording');%%This command is finish recording
% play(a);%Play recorded audio.
% b=getaudiodata(a);%Store recorded audio signal in numeric array
% audiowrite('sss11.wav',b,Fs); % write the sound file wav. format
% ii.
[y,Fs]=audioread('C:\Users\TOSHIBA\Desktop\dspproject\sss11.wav'); %this line reads audio signal from its location on the computer
info = audioinfo('C:\Users\TOSHIBA\Desktop\dspproject\sss11.wav') %this command gives information about the sound being recorded
% iii.
Ts=1/Fs;%define sampling period
t = 0:Ts:(length(y)-1)*Ts; % define time vector
fft_y = fft(y); %this line receives the FFT of the original signal
f=linspace(0,Fs,length(fft_y));% define frequcy vector length
n =floor(length(f)/2); %this line got half the length of the vector and we rounded it up to the smallest number
fk=f/1000;%this line convert to frequency kHz unit
c=mag2db(abs(fft_y));%this line convert the magnitude to a logarithmic value(dB).
figure;%create new figure window (for figure1)
subplot(211);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),c(1:n));%this line plot the magnitude response of voice data
title('Magnitude response of voice data');%write a title for magnitude response of voice data graph
xlabel('frequency(kHz)');%write a x axis name for magnitude response of voice data graph
ylabel('magnitude(dB)');%write a y axis name for magnitude response of voice data graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
phase1=angle((fft_y));%the angle of the original signal that received the fast fourier transform is defined
phi1=unwrap(phase1);%with the unwrap command, this angle is made continuous
subplot(212);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),phi1(1:n));%this line plot the phase response of voice data
title('Phase response of voice data');%write a title for phase response of voice data graph
xlabel('frequency(kHz)');%write a x axis name for phase response of voice data graph
ylabel('phase(deg)');%write a y axis name for phase response of voice data graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% iv.
ty=transpose(y);%the transpose of the y vector is taken from
fc=10e3;  %define 10kHz frequency carrier signal
w=2*pi*fc;%define w value
s=2*max(y)*sin(w*t); %create tone signal
 % v
d=s+ty; %Add tone signal to voice signal and obtain distorted signal
% vi.
fft_d = fft(d); %this line receives the FFT of the distorted signal.
d1=mag2db(abs(fft_d));%convert to distorted signal magnitude to decibel unit
figure;%create new figure window (for figure2)
subplot(211);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),d1(1:n));%this line plot the magnitude response of distorted signal
title('Magnitude response of distorted signal');%write a title for magnitude response of distorted signal graph
xlabel('frequecy (kHz)');%write a x axis name for magnitude response of distorted signal graph
ylabel('magnitude(dB)');%write a y axis name for magnitude response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
phase2=angle(fft_d);%the angle of the distorted signal that received the fast fourier transform is defined
phi2=unwrap(phase2);%with the unwrap command, this angle is made continuous
subplot(212);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),phi2(1:n));%this line plot the phase response of distorted signal
title('Phase response of distorted signal');%write a title for phase response of distorted signal graph
xlabel('frequency(kHz)');%write a x axis name for phase  response of distorted signal graph
ylabel('phase(deg)');%write a y axis name for phase response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% vii.
figure;%create new figure window (for figure3)
subplot(311);%in the figure window, it divides into rectangles, and 3 is Row, 1 is column, and 1 is the rectangle to be placed in.
plot(t,y);%this line plot the original signal
title('original signal');%write a title for original signal graph
xlabel('time(sec)');%write  x axis name for original signal graph
ylabel('amplitude(v)');%write  y axis name for original signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(312);%in the figure window, it divides into rectangles, and 3 is Row, 1 is column, and 2 is the rectangle to be placed in.
plot(t,d);%this line plot the distorted signal
title('distorted signal');%write a title for distorted signal graph
xlabel('time(sec)');%write  x axis name for distorted signal graph
ylabel('amplitude(v)');%write  y axis name for distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(313);%in the figure window, it divides into rectangles, and 3 is Row, 1 is column, and 3 is the rectangle to be placed in.
plot(t,s);%this line plot the tone signal
title('tone signal');%write a title for tone signal graph
xlabel('time(sec)');%write  x axis name for tone signal graph
ylabel('amplitude(v)');%write  y axis name for tone signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% viii.
figure;%create new figure window (for figure4)
subplot(221);%in the figure window, it divides into rectangles, and 2 is Row, 2 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),c(1:n)); %this line plot the magnitude response of original signal
title('Magnitude response of original signal'); %write a title for magnitude response of orginal signal graph
xlabel('frequency(kHz)'); %write  x axis name for magnitude  response of original signal graph
ylabel('magnitude(dB)'); %write a y axis name for magnitude  response of original signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(222);%in the figure window, it divides into rectangles, and 2 is Row, 2 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),phi1(1:n));%this line plot the phase response of original signal
title('Phase response of original signal');%write a title for phase response of orginal signal graph
xlabel('frequency(kHz)');%write  x axis name for phase  response of original signal graph
ylabel('phase(deg)'); %write a y axis name for phase  response of original signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(223);%in the figure window, it divides into rectangles, and 2 is Row, 2 is column, and 3 is the rectangle to be placed in.
plot(fk(1:n),d1(1:n));%this line plot the magnitude response of distorted signal
title('Magnitude response of distorted signal');%write a title for magnitude response of distorted signal graph
xlabel('frequecy (kHz)');%write  x axis name for magnitude response of distorted signal graph
ylabel('magnitude(dB)');%write  y axis name for magnitude response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(224);%in the figure window, it divides into rectangles, and 2 is Row, 2 is column, and 4 is the rectangle to be placed in.
plot(fk(1:n),phi2(1:n));%this line plot the phase response of distorted signal
title('Phase response of distorted signal');%write a title for phase response of distorted signal graph
xlabel('frequency(kHz)');%write  x axis name for phase response of distorted signal graph
ylabel('phase(deg)');%write  y axis name for phase response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% ix.
% Designing minimum order FIR-LPF with PassbandFrequency,StopbandFrequency,PassbandRipple and StopbandAttenuation
%  And Specifying Frequency Parameters in Hertz
 df = designfilt('lowpassfir','PassbandFrequency',9700,...
   'StopbandFrequency',10000,'PassbandRipple',0.17,...
   'StopbandAttenuation',70,'SampleRate',44100);
% x. 
%(for figure5)
 hfvt = fvtool(df,'Analysis','freq');% this line overlays the magnitude and phase response of the designed filter.
 % xi.
x = filter(df,d);% This step inserts the distorted signal into the designed filter, resulting in the X parameter at its output.
figure;%create new figure window (for figure6)
plot(t,x);%this line plot the filtered signal
title('filtered signal');%write a title for filtered signal graph
xlabel('time(sec)');%write  x axis name for filtered signal graph
ylabel('amplitude(v)');%write  y axis name for filtered signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xii.
fft_x = fft(x); %this line receives the FFT of the filtered signal.
x1=mag2db(abs(fft_x));%this line converted to filtered signal magnitude to decibel unit
figure;%create new figure window (for figure7)
subplot(211);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),x1(1:n));%this line plot the magnitude response of distorted signal
title('Magnitude response of filtered signal');%write a title for magnitude response of filtered signal graph
xlabel('frequency(kHz)');%write  x axis name for magnitude response of filtered signal graph
ylabel('Magnitude(dB)');%write  y axis name for magnitude response of filtered signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
phasex=angle(fft_x);%the angle of the filtered signal that received the fast fourier transform is defined
phix=unwrap(phasex);%with the unwrap command, this angle is made continuous
subplot(212);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),phix(1:n));%this line plot the phase response of distorted signal
title('Phase response of filtered signal');%write a title for phase response of filtered signal graph
xlabel('frequency(kHz)');%write  x axis name for phase response of filtered signal graph
ylabel('Phase(deg)');%write  y axis name for phase response of filtered signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xiii.
subplot(231);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),c(1:n));%this line plot the magnitude response of voice data
title('Magnitude response of voice data');%write a title to magnitude response of voice data graph
xlabel('frequency(kHz)');%write  x axis name to magnitude response of voice data graph
ylabel('magnitude(dB)');%write  y axis name to magnitude response of voice data graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(234);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 4 is the rectangle to be placed in.
plot(fk(1:n),phi1(1:n));%this line plot the phase response of voice data
title('Phase response of voice data');%write a title to phase response of voice data graph
xlabel('frequency(kHz)');%write  x axis name to phase response of voice data graph
ylabel('phase(deg)');%write  y axis name to phase response of voice data graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(232);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),d1(1:n));%this line plot the magnitude response of distorted signal
title('Magnitude response of distorted signal');%write a title to magnitude response of distorted signal
xlabel('frequecy (kHz)');%write  x axis name to magnitude response of distorted signal graph
ylabel('magnitude(dB)');%write  y axis name to magnitude response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(235);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 5 is the rectangle to be placed in.
plot(fk(1:n),phi2(1:n));%this line plot the phase response of distorted signal
title('Phase response of distorted signal');%write a title to phase response of distorted signal
xlabel('frequency(kHz)');%write  x axis name to phase response of distorted signal graph
ylabel('phase(deg)');%write  y axis name to phase response of distorted signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(233);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 3 is the rectangle to be placed in.
plot(fk(1:n),x1(1:n));%this line plot the magnitude response of filtered signal
title('Magnitude response of filtered signal');%write a title for magnitude response of filtered signal
xlabel('frequency(kHz)');%write  x axis name for magnitude response of filtered  signal graph
ylabel('Magnitude(dB)');%write  y axis name for magnitude response of filtered  signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(236);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 6 is the rectangle to be placed in.
plot(fk(1:n),phix(1:n));%this line plot the phase response of filtered signal
title('Phase response of filtered signal');%write a title for phase response of filtered signal
xlabel('frequency(kHz)');%write  x axis name for phase response of filtered  signal graph
ylabel('Phase(deg)');%write  y axis name for phase response of filtered  signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xiv.
figure;%create new figure window (for figure8)
plot(t,y);%this line plot the original signal
hold on;% this command draws two different signals on the same axis
xlim([0.90,0.95]);%this command limits the time axis to 0.05 seconds.
plot(t,x);%this line plot the filtered signal
title('original signal and filtered signal');%write a title for these graphs
xlabel('time(sec)');%write a x axis name for these graphs
ylabel('amplitude(v)');%write a y axis name for these graphs
legend('original-s','filtered-s');% this command used to distinguish graphs drawn on the same axis
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xv.
% filter delays some frequency components more than others.(phase distortion).
% because of the phase distortion delay can occur.
D = mean(grpdelay(df))%writing delay at command window.
xs = filter(df,[y; zeros(D,1)]); % Append D zeros to the input data
xs = xs(D+1:end);         % Shift data to compensate for delay
% xvi.
fft_xs = fft(xs); %this line receives the FFT of the time delay compensated signal.
xs1=mag2db(abs(fft_xs));%this line  converted the magnitude to a logarithmic value.
figure;%create new figure window (for figure9)
subplot(211);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),xs1(1:n));%this line plot the magnitude response of time delay compensated signal 
title('Magnitude response of time delay compensated signal');%write a title for magnitude response of time delay compensated signal graph
xlabel('frequency(kHz)');%write a x axis name for the magnitude response of time delay compensated signal graph
ylabel('magnitude(dB)');%write a y axis name for the magnitude response of time delay compensated signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
phases=angle((fft_xs));%the angle of the time delay compensated signal that received the fast fourier transform is defined
phis=unwrap(phases);%with the unwrap command, this angle is made continuous
subplot(212);%in the figure window, it divides into rectangles, and 2 is Row, 1 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),phis(1:n));%this line plot the phase response of time delay compensated signal 
title('Phase response of time delay compensated signal');%write a title for phase response of time delay compensated signal graph
xlabel('frequency(kHz)');%write a x axis name for the phase response of time delay compensated signal graph
ylabel('phase(deg)');%write a y axis name for the phase response of time delay compensated signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xvii.
figure;%create new figure window (for figure10)
subplot(231);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 1 is the rectangle to be placed in.
plot(fk(1:n),c(1:n)); %this line plot the magnitude response of original signal
title('Magnitude response of original signal'); %write a title for magnitude response of original signal graph
xlabel('frequency(kHz)'); %write  x label name for magnitude response of original signal graph
ylabel('magnitude(dB)'); %write y label name for magnitude response of original signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(234);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 4 is the rectangle to be placed in.
plot(fk(1:n),phi1(1:n));%this line plot the phase response of original signal
title('Phase response of original signal'); %write a title for phase response of original signal graph
xlabel('frequency(kHz)');%write  x label name for phase response of original signal graph
ylabel('phase(deg)');%write y label name for phase response of original signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(232);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 2 is the rectangle to be placed in.
plot(fk(1:n),x1(1:n));%this line plot the magnitude response of filtered signal
title('Magnitude response of filtered signal');%write a title for magnitude response of filtered signal graph
xlabel('frequency(kHz)');%write x label name for magnitude response of filtered signal graph
ylabel('Magnitude(dB)');%write y label name for magnitude response of filtered signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(235);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 5 is the rectangle to be placed in.
plot(fk(1:n),phix(1:n));%this line plot the phase response of filtered signal
title('Phase response of filtered signal');%write a title for phase response of filtered signal graph
xlabel('frequency(kHz)');%write x label name for phase response of filtered signal graph
ylabel('Phase(deg)');%write y label name for phase response of filtered signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(233);%in the figure window, it divides into rectangles, and 2 is Row, 3 is column, and 6 is the rectangle to be placed in.
plot(fk(1:n),xs1(1:n));%this line plot the magnitude response of time delay compensated signal
title('Magnitude response of time delay compensated signal');%write a title for magnitude response of time delay compensated signal graph
xlabel('frequency(kHz)');%write x label name for magnitude response of  time delay compensated signal graph
ylabel('magnitude(dB)');%write y label name for magnitude response of  time delay compensated signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
subplot(236);%in the figure window, it divides into rectangles, and 2 is Row, 2 is column, and 3 is the rectangle to be placed in.
plot(fk(1:n),phis(1:n));%this line plot the phase response of time delay compensated signal
title('Phase response of time delay compensated signal');%write a title for phase response of time delay compensated signal graph
xlabel('frequency(kHz)');%write x label name for phase response of  time delay compensated signal graph
ylabel('phase(deg)');%write y label name for phase response of  time delay compensated signal graph
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% xviii.
figure; %create new figure window (for figure11)
subplot(211);
plot(t,y);%this line plot the original signal
hold on;% this command draws two different signals on the same axis
xlim([0.90,0.95]);%this command limits the time axis to 0.05 seconds.
plot(t,x);%this line plot the filtered signal
title('original signal and filtered signal');%write a title for these graphs
xlabel('time(sec)');%write a x axis name for these graphs
ylabel('amplitude(v)');%write a y axis name for these graphs
legend('original-s','filtered-s');% this command used to distinguish graphs drawn on the same axis
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
% 
subplot(212);
plot(t,y);%this line plot the original signal
hold on;% this command draws two different signals on the same axis
xlim([0.90,0.95]);%this command limits the time axis to 0.05 seconds.
plot(t,xs);%this line plot the time delay compensatd filtered signal
title('original signal and timedelay-com-filtered signal');%write a title for these graphs
xlabel('time(sec)');%write a x axis name for these graphs 
ylabel('amplitude(v)');%write a y axis name for these graphs
legend('original-s','timedelay-compansated-filtered-s');%this command used to distinguish graphs drawn on the same axis
grid on;% this command divides the figure windows into squares, which increases the readability of the graph.
 % xix.
% sound(y,Fs)  %play original signal
% pause(5);%stop playing original signal
% sound(d,Fs)  %play distorted signal
% pause(5);%stop playing distorted signal
% sound(x,Fs)  %play filtered signal
% pause(5);%stop playing filtered signal
% sound(xs,Fs) %play time delay compensated filtered signal
% pause(5);%stop playing time delay compensated filtered signal
% sound(s,Fs)  %play 10kHz tone signal
% pause(5);%stop playing 10kHz tone signal
