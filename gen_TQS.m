function img_seq = gen_TQS(n_display_per_image, n_rsvp)
% generate a trial queue sequence
% n_display_per_image is a vector with the length of img_fnames and each
% element corresponds number of requested presentations per image

n_trs_tot = ceil(sum(n_display_per_image)/n_rsvp); 
if mod(n_trs_tot,n_rsvp) ~=0
    n_trs_tot = n_trs_tot + (n_rsvp - mod(n_trs_tot,n_rsvp));
end
seq = repelem(linspace(1,length(n_display_per_image),length(n_display_per_image)), n_display_per_image); 

seq = seq(randperm(length(seq))); 
if length(seq) < n_trs_tot * n_rsvp
    seq = [seq, randperm(length(n_display_per_image), n_trs_tot*n_rsvp - length(seq))];
end

img_seq = reshape(seq,[n_rsvp,n_trs_tot]); 
end
