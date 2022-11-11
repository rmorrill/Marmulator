function read_lick_cmd = gen_lickometer_command(lick_arduino)


lick_trig_pin = lick_arduino.lick_pin; 
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
read_lick_cmd = ['1' Alphabet(lick_trig_pin + 1)];


