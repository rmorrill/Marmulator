function pix_sz = convertDegreesPix(obj_deg, screen_dist, screen_sz, screen_pix)
deg_per_pix = convertPixDegrees(1, screen_dist, screen_sz, screen_pix); 
pix_per_deg = 1/deg_per_pix; 
pix_sz = obj_deg * pix_per_deg; 
