function assign_pump_pins(pump_arduino)

ahand = pump_arduino.ahand;
pump_trig_pin = pump_arduino.pump_pin; 
pump_led_trig_pin = pump_arduino.pump_led_pin; 

IOPort('Flush', ahand);
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
pump_trig_letter = Alphabet(pump_trig_pin + 1);
pump_led_trig_letter = Alphabet(pump_led_trig_pin + 1);

assign_pumppin_cmd = ['0' pump_trig_letter '1']';
assign_pumpledpin_cmd = ['0' pump_led_trig_letter '1']'; 

fprintf('assign pump cmd: %s\n', assign_pumppin_cmd)
fprintf('assign pump led cmd: %s\n', assign_pumpledpin_cmd)


IOPort('Write', ahand, assign_pumppin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_pumppin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to pump trigger\n', pump_trig_pin); 

IOPort('Write', ahand, assign_pumpledpin_cmd, 1);
WaitSecs(0.05);
IOPort('Write', ahand, assign_pumpledpin_cmd, 1);
WaitSecs(0.05); 
fprintf('Assigned pin %d to pump led trigger\n', pump_led_trig_pin); 


