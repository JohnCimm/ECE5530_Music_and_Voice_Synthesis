%% Music Synthesis Lab 4 -- Ian Ashworth, Emily Erickson, John Cimmarusti


%% 4.a Determine Sampling Frequency
fs = 8000;

%% 4.b Determine note duration and frequency

load bach_fugue.mat
num_melodies = length(theVoices); % Tells the number of melodies present
% Creates vectors containing the keynumber for each note
melody_1_notes = theVoices(1).noteNumbers; 
melody_2_notes = theVoices(2).noteNumbers;
melody_3_notes = theVoices(3).noteNumbers;
% Determines the time in seconds a beat is
bpm = 120;
bps = bpm/60;
spb = 1/bps;
spp = spb/4;

% Gets the difference between each note and note 49 which has a freq of 440
% Hz
for i = 1:length(melody_1_notes)
    m1_freq_ratio(i) = 49-melody_1_notes(i);
end
for i = 1:length(melody_2_notes)
    m2_freq_ratio(i) = 49-melody_2_notes(i);
end
for i = 1:length(melody_3_notes)
    m3_freq_ratio(i) = 49-melody_3_notes(i);
end

% Using the ratios the freqency of each note is calculated and stored
for i = 1:length(melody_1_notes)
    if m1_freq_ratio(i) < 0
        m1_note_freq(i) = 440/2^(abs(m1_freq_ratio(i))/12);
    elseif m1_freq_ratio(i) > 0
         m1_note_freq(i) = 440*2^(abs(m1_freq_ratio(i))/12);
    else
        m1_note_freq(i) = 440;
    end
end

m1_note_freq

for i = 1:length(melody_2_notes)
    if m2_freq_ratio(i) < 0
        m2_note_freq(i) = 440/2^(abs(m2_freq_ratio(i))/12);
    elseif m2_freq_ratio(i) > 0
         m2_note_freq(i) = 440*2^(abs(m2_freq_ratio(i))/12);
    else
        m2_note_freq(i) = 440;
    end
end

for i = 1:length(melody_3_notes)
    if m3_freq_ratio(i) < 0
        m3_note_freq(i) = 440/2^(abs(m3_freq_ratio(i))/12);
    elseif m3_freq_ratio(i) > 0
         m3_note_freq(i) = 440*2^(abs(m3_freq_ratio(i))/12);
    else
        m3_note_freq(i) = 440;
    end
end

% Calculate and stores each note duration in the different melodies

melody_1_note_duration = zeros(size(melody_1_notes));
for i = 1:length(melody_1_notes)
    melody_1_note_duration(i) = theVoices(1).durations(i)*spp;
end

melody_1_note_duration

melody_2_note_duration = zeros(size(melody_2_notes));
for i = 1:length(melody_2_notes)
    melody_2_note_duration(i) = theVoices(2).durations(i)*spp;
end
    
    melody_3_note_duration = zeros(size(melody_3_notes));
for i = 1:length(melody_3_notes)
    melody_3_note_duration(i) = theVoices(3).durations(i)*spp;
end

%% 4.c Music synthesis

m1synth = zeros(1, round(sum(melody_1_note_duration)*fs+length(melody_1_notes)) );
n1 = 1;
for kk = 1:length(melody_1_notes)
    tt = 0:1/fs:melody_1_note_duration(kk);
    freq = m1_note_freq(kk);
    tone = real(1.*exp(j*2*pi*freq*tt));
    n2 = n1 + length(tone)-1;
    Ltone = length(tone);
    tail = round( max([0.04*fs,0.05*Ltone]) );
    m1synth(n1:n2) = m1synth(n1:n2) +  mg_env(Ltone,tail).*tone;
    n1 = n2+1;
end
audiowrite('First_Melody_Synth.wav',m1synth,fs)
soundsc(m1synth,fs);
%% 4.d Plotting sinusoids

figure
for i = 1:3
    tt = 0:1/fs:(melody_1_note_duration(i));
    freq = m1_note_freq(i);
    tone = real(1.*exp(j*2*pi*freq*tt));
    Ltone = length(tone);
    tail = round( max([0.04*fs,0.05*Ltone]) );
    subplot(3,1,i)
    plot(tt,mg_env(Ltone,tail).*tone)
    title("Note "+ i)
    xlabel("Time (s)")
    ylabel("Amplitude")
    xline([0:1/freq:melody_1_note_duration(i)])
end

% The black lines on the graph represent each period for the note freqeuncy
% for the first three notes.  Notes 1 and 3 are the same, but note 3 has a
% longer duration so it has more periods.

%% 4.e Spectrograph

figure
subplot(2,1,1)
specgram(m1synth,512,fs)
title("Synthesized First Melody")
ylim([0,1400])
subplot(2,1,2)
plot(m1_note_freq)
xlim([0,123])
xlabel("Melody Note Number")
ylabel("Note Frequency (Hz)")
title("Calculated Note Frequency for First Melody")



%% Other Melodies

m2synth = zeros(1, round(sum(melody_2_note_duration)*fs+length(melody_2_notes)) );
n1 = 1;
for kk = 1:length(melody_2_notes)
    tt = 0:1/fs:melody_2_note_duration(kk);
    freq = m2_note_freq(kk);
    tone = real(1.*exp(j*2*pi*freq*tt));
    n2 = n1 + length(tone)-1;
    Ltone = length(tone);
    tail = round( max([0.04*fs,0.05*Ltone]) );
    m2synth(n1:n2) = m2synth(n1:n2) +  mg_env(Ltone,tail).*tone;
    n1 = n2+1;
end


m3synth = zeros(1, round(sum(melody_3_note_duration)*fs+length(melody_3_notes)) );
n1 = 1;
for kk = 1:length(melody_3_notes)
    tt = 0:1/fs:melody_3_note_duration(kk);
    freq = m3_note_freq(kk);
    tone = real(1.*exp(j*2*pi*freq*tt));
    n2 = n1 + length(tone)-1;
    Ltone = length(tone);
    tail = round( max([0.04*fs,0.05*Ltone]) );
    m3synth(n1:n2) = m3synth(n1:n2) +  mg_env(Ltone,tail).*tone;
    n1 = n2+1;
end


function envelope = mg_env(length, attack)
    envelope = ones(1,length);
    attack1 = max(2,round(attack/2));
    envelope(1:attack1) = (0:(attack1-1)/attack1);
    release = attack;
    e(length:-1:length-release) = (0:release)/release;
end

