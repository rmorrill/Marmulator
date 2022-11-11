function assign_lickometer_pins(lick_arduino)
% three types of triggers: 
% trial: goes on at beginning of trial (inlcudes pre-stim initiation period) 
% stimulus: goes on at presentation of stimulus, including for each RSVP
% session: goes on at beginning of behavior session and off after all
% trials or exist

ahand = lick_arduino.ahand;
lick_trig_pin = lick_arduino.lick_pin; 

IOPort('Flush', ahand);
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
lick_trig_letter = Alphabet(lick_trig_pin + 1);

assign_lickpin_cmd = ['0' lick_trig_letter '0']';


IOPort('Write', ahand, assign_lickpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_lickpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to lickometer trigger\n', lick_trig_pin); 


