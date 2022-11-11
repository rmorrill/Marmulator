function checkLickometer_gui(obj, event, ahand, read_lick_cmd, lick_box_hand)

IOPort('Write', ahand, read_lick_cmd, 2);
val = IOPort('Read', ahand, 1,3);
lick = val(1)-48;
lick_col = [0.07, 0.62, 1.00]; 
orig_col = [0.5, 0.5, 0.5]; 
if lick
    %disp('Lick'); 
    set(lick_box_hand, 'BackgroundColor', lick_col); 
    drawnow; 
    WaitSecs(0.05);
    set(lick_box_hand, 'BackgroundColor', orig_col); 
    drawnow; 
end

