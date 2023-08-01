function write_pump_cmd = gen_pump_command(pump_arduino)


pump_trig_pin = pump_arduino.pump_pin; 
Alphabet = 'abcdefghijklmnopqrstuvwxyz';
pump_pin_letter= Alphabet(pump_trig_pin + 1);


write_pump_cmd.on = ['2' pump_pin_letter '1'];
write_pump_cmd.off = ['2' pump_pin_letter '0'];


