#!/usr/bin/perl

for($rmm=0;$rmm<4;$rmm++){
    for($tfb=0;$tfb<48;$tfb++){
	for($trip=0;$trip<4;$trip++){
	    for($ch=0;$ch<16;$ch++){
		system "./dump_INGRID_Ch_config $rmm $tfb $trip $ch";
	    }
	}
    }
}
