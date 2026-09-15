function assign_2p_trig_pins(trig_arduino)

ahand = trig_arduino.ahand;
trig_2p_pin = trig_arduino.trig_2p_read_pin; 

IOPort('Flush', ahand);
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
trig_2p_letter = Alphabet(trig_2p_pin + 1);

assign_trigpin_cmd = ['0' trig_2p_letter '0']';


IOPort('Write', ahand, assign_trigpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_trigpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to 2p trigger\n', trig_2p_pin); 