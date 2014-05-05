#!/usr/bin/perl -w
use strict;


# Usage $ program.pl <uniq_snps> <data_file>
# uniq_snps:
# rs0001
# rs0002
# rs0006
# rs0001
# ...


my $uniq_filename = shift(@ARGV);
open(my $fh, "<$uniq_filename");

my %snps;
while(<$fh>)
{
    chomp;
    my $line = $_;
    # $line =~ s/\[\[\:space\:\]\]//;
    $snps{$line} += 1;
}
close($fh);

while(<>)
{
    # rs0123   USA  A   C   1\n
    chomp;
    # rs0123   USA  A   C   1
    my $line = $_;
    # rs0123   USA  A   C   1
    my @split_line = split(/\s+/, $line);
    # array = ["rs0123", "USA", "A", "C", "1"]
    if (defined($snps{$split_line[0]}) || ($split_line[0]=~/[SNP]/))
    {
        print "$line\n";
    }
}

