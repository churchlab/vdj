#!/usr/bin/perl

open(FNAFILE, $ARGV[0]);
open(QUALFILE, $ARGV[1]);

while($fheader = <FNAFILE>){
    $qheader = <QUALFILE>;
    chomp($fheader);
    chomp($qheader);

    if($fheader == $qheader){
        $fna = <FNAFILE>;
        $qual = <QUALFILE>;
        chomp($fna);
        chomp($qual);
        $fheader =~ s/>/@/;
        print STDOUT $fheader . "\n" . $fna . "\n+\n" . $qual . "\n";
    }

}
