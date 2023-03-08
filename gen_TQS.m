function img_seq = gen_TQS(n_display_per_image, n_rsvp)
% generate a trial queue sequence

n_trs_tot = ceil(sum(n_display_per_image)/n_rsvp); 

seq = repelem(linspace(1,length(n_display_per_image),length(n_display_per_image)), n_display_per_image); 

seq = seq(randperm(length(seq))); 
if length(seq) < n_trs_tot
    seq = [seq, randperm(length(n_display_per_image), n_trs_tot - length(seq))];
end

img_seq = reshape(seq,[n_rsvp,n_trs_tot]); 
end
