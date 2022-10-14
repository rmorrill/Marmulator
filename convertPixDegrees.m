function deg_sz = convertPixDegrees(obj_pix, screen_dist, screen_sz, screen_pix)

pix_per_cm = screen_pix(1)/screen_sz(1); 
obj_cm = obj_pix/pix_per_cm; 
deg_sz = rad2deg(2*atan((obj_cm/2)/screen_dist)); 
