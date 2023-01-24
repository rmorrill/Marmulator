function assign_trigger_pins(trigger_arduino)
% three types of triggers: 
% trial: goes on at beginning of trial (inlcudes pre-stim initiation period) 
% stimulus: goes on at presentation of stimulus, including for each RSVP
% session: goes on at beginning of behavior session and off after all
% trials or exist
% sample command: per elias's request, goes on when the monkey initiates a
% trial and off after a reward/punish screen (2023-01-19 YJ) 

ahand = trigger_arduino.ahand;
trial_trig_pin = trigger_arduino.trial_pin;
stim_trig_pin = trigger_arduino.stim_pin;
session_trig_pin = trigger_arduino.session_pin;
sc_trig_pin = trigger_arduino.sc_pin; 

trigger_flag = true;
IOPort('Flush', ahand);
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
trial_pin_letter = Alphabet(trial_trig_pin + 1);
stim_trig_letter = Alphabet(stim_trig_pin + 1);
session_trig_letter = Alphabet(session_trig_pin + 1);
sc_trig_letter = Alphabet(sc_trig_pin + 1); 

assign_trialpin_cmd = ['0' trial_pin_letter '1']';
assign_stimpin_cmd = ['0' stim_trig_letter '1']';
assign_sesspin_cmd = ['0' session_trig_letter '1']';
assign_scpin_cmd = ['0' sc_trig_letter '1']; 

IOPort('Write', ahand, assign_sesspin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_sesspin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to session trigger\n', session_trig_pin); 

IOPort('Write', ahand, assign_trialpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_trialpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to trial trigger\n', trial_trig_pin); 

IOPort('Write', ahand, assign_stimpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_stimpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to stimulus trigger\n', stim_trig_pin); 

IOPort('Write', ahand, assign_scpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_scpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to sample command trigger\n', sc_trig_pin); 



