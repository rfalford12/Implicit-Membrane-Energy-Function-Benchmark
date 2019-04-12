#!/usr/bin/perl

$pdb = $ARGV[0];
$pdb = lc($pdb);
$id = substr($pdb, 1, 2);

system ("curl \'http://pdbtm.enzim.hu/data/database/".$id."/".$pdb.".trpdb.gz\' --output ".uc($pdb)."_tr.pdb.gz");
system ("gunzip ".uc($pdb)."_tr.pdb.gz");

