#!/usr/bin/perl
for($noise=2; $noise<20; $noise=$noise+2){
    system "./IngAddNoiseMC -r 111 -s 3 -m -n $noise";
    system "./IngTClster -f ~/data/mc_data/add_noise/ingmc_numu_fe_nd3_horn250ka_111_noiseadd${noise}.root";
    system "./IngBasicRecon -r 77 ";
    system "./IngReconMuBasic -r 77 ";
    system "mv ~/data/dst/ingrid_00000077_0000_track.root ~/data/mc_data/add_noise/mc_afteraddnoise_${noise}.root";
}
