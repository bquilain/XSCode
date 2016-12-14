#!/bin/zsh +x
for i in {1..1000};do
#    /home/bquilain/CC0pi_XS/XS/UnfoldingOptimisation_Dvt -d /home/bquilain/CC0pi_XS/XS/files/MCSelected_Systematics0_0 -m /home/bquilain/CC0pi_XS/XS/files/MCSelected_Systematics0_0 -o /home/bquilain/CC0pi_XS/XS/files/MCUnfolded_Statistics${i} -n 1 -s 1
    mv /home/bquilain/CC0pi_XS/XS/files/MCUnfolded_Statistics${i}.root /home/bquilain/CC0pi_XS/XS/files/MCUnfolded_Systematics1_${i}.root
done
