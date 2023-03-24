%% Main script

im_tag = 'Washington_4'; % im_tag = 'Fields';
demo_formation_mrca(im_tag);
demo_formation_classic(im_tag);
demo_formation_cassi(im_tag);

im_tag = 'Janeiro'; % im_tag = 'Stockholm';
demo_reconstruction_jodefu(im_tag);
demo_reconstruction_classic(im_tag);

im_tag = 'Beijing_4';
demo_parameters(im_tag);