function [sess_trig_cmd, trial_trig_cmd, stim_trig_cmd, sc_trig_cmd] = gen_trig_commands(trigger_arduino)

trial_trig_pin = trigger_arduino.trial_pin;
stim_trig_pin = trigger_arduino.stim_pin;
session_trig_pin = trigger_arduino.session_pin;
sc_trig_pin = trigger_arduino.sc_pin; 

Alphabet = 'abcdefghijklmnopqrstuvwxyz';
trial_pin_letter = Alphabet(trial_trig_pin + 1);
stim_trig_letter = Alphabet(stim_trig_pin + 1);
session_trig_letter = Alphabet(session_trig_pin + 1);
sc_trig_letter = Alphabet(sc_trig_pin + 1); 

sess_trig_cmd.on = ['2' session_trig_letter '1'];
sess_trig_cmd.off = ['2' session_trig_letter '0'];

trial_trig_cmd.on = ['2' trial_pin_letter '1'];
trial_trig_cmd.off = ['2' trial_pin_letter '0'];

stim_trig_cmd.on = ['2' stim_trig_letter '1'];
stim_trig_cmd.off = ['2' stim_trig_letter '0'];

sc_trig_cmd.on = ['2' sc_trig_letter '1'];
sc_trig_cmd.off = ['2' sc_trig_letter '0']; 

