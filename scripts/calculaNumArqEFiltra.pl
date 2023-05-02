#!/usr/bin/perl
use strict;
my $bin_size = shift(@ARGV);
#my $locus_bin_size = shift(@ARGV);
my $threshold = shift(@ARGV); # we establish a limit on the number or architectures
my %histogram_number;
my %histogram_locus;
my $max = 0;
my $max_locus = 0;
while (my $line = <>){
  chomp($line);
  (my $name, my $tfbs) = split(/\t/,$line);
  print STDERR "name:$name,tfbs=$tfbs\n";
  #my @locuses = split(",", $tfbs);
  my @locuses = split(";", $tfbs);
  my $number = 2**(scalar(@locuses));
  #foreach my $locus (@locuses){
  #  my @entries = split(/\//,$locus);
  #  $number *= scalar(@entries);
  #}
  if ($number > $max){$max = $number};
  #my $tfbs_num = scalar(split(/\,|\//,$tfbs));
  #if ($tfbs_num > $max_locus){$max_locus=$tfbs_num};
  my $bin = (int($number/$bin_size)+1) * $bin_size;
  $histogram_number{$bin} ++;
  #print STDERR "$name,$number\n";
  #$bin = (int($tfbs_num/$locus_bin_size)+1) * $locus_bin_size;
  #$histogram_locus{$bin} ++;
  #print STDERR "$name,$number\n";
  if ($number <= $threshold){ 
    print $line,"\n";
  }
  else{
    print STDERR "WARNING: Gene with too many subsets ($number):\n\t$line\n";
  }
}
print STDERR "Subset Number Histogram:\n"; 
foreach my $key (sort {$a <=> $b} (keys(%histogram_number))){
  print STDERR "$key,$histogram_number{$key}\n";
}
#print STDERR "Locus Number Histogram:\n"; 
#foreach my $key (sort {$a <=> $b} (keys(%histogram_locus))){
#  print STDERR "$key,$histogram_locus{$key}\n";
#}
print "max=$max,max ";#locus=$max_locus\n"; 
