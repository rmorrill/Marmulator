function timerErrorFcn_gui(obj, event)

t = timerfind('Name', 'lickometer_gui_timer'); 

fprintf('******** TIMER ERROR! lickometer_gui_timer has been stopped.\nDid you "clear all"? Use "clearvars" instead. To restart, click lickometer Disconnect and then Connect\n'); 
