#!/usr/local/bin/perl
use Getopt::Long;

# to call this program use "perl subsets2.pl[-t treefile] [-r] [-l]
use strict;
use Set::Scalar;
use Set::Bag;
#
#we will assume TFBS description includes a name, followed by start and end positions, as in "AB8(1234:1254)
#

#
#WARNING, FOR PATTERN REPORT WE IGNORE POSITIONS BASED ON RESULTS OF
#
#
my %patternHash;
my %arch_hash;
my %geneGroupHash;
my %geneFamilyHash;
GetOptions( 'treefile=s' => \my $treefile #if we have argument -t loads treefile parameter
          , 'repeats' => \my $repeats  # -r and --repeats indicate we consider repeats
          , 'locus' => \my $locus  # -l and --locus indicate we want to uses locuses instead of individual tfbs
          );
#first open tree file, and parse tree, building a hash with
#all gene groups
my $expression_cluster_hash;
if ($treefile){
    open(TREE,$treefile) or die "could not open expression clustering file $treefile\n";
    my $tree_text = <TREE>;
    $expression_cluster_hash = buildSetsFromTree($tree_text);
}
#now we will compute the reports


my $pattern_count;
while (my $line = <>){
  chomp($line);
  print STDERR "reading=>$line<=\n";
  (my $gene_name, my $tfbs_string) = split(/\t/,$line);
  #$tfbs_string =~ s/\,/;/g;
  $patternHash{$tfbs_string} .= "$gene_name\n";
  my $subset_string = subsets($tfbs_string, $locus,-1);#obs: -1 indicates no pre-restrictions on the position
  print STDERR "subsets:$subset_string\m"; 
  #Next step cleans the architectures removing repeats and order
  my $clean_subsets = clean($subset_string,$repeats, $locus);
  print STDERR "dirty_subset:$subset_string;\nclean subsets:$clean_subsets\n";
  #my @architectures = split(";", $subset_string);
  my @architectures = split(";", $clean_subsets);
  #now, for each architecture found for the gene, add the
  #gene in the architecture entry of $arch_hash
  foreach my $one_arch (@architectures){
      print  STDERR "\tadding architecture:$one_arch\n";
      if ($arch_hash{$one_arch} !~ m/$gene_name/){
	  $arch_hash{$one_arch} .= "$gene_name,";
      }
      else {
	#this should not happen, if it does, there is a problem in subsets() or clean()
	print "\tWARNING one architecture twice in a gene\n";
      }
    }
  #checking the number of architectures for this gene.
  $pattern_count += scalar(split(";", $subset_string));
  print STDERR "\tnumber of subsets", scalar(split(";", $clean_subsets)),"\n";
  print STDERR "\tsubsets:$clean_subsets\n";
}
#will print hash in decreasing order of number of genes for a pattern
print "ARCHITECTURE REPORT:\n==================\n\n";
foreach my $key (sort {my $numa = scalar(split(",",$arch_hash{$a}));
		       my $numb = scalar(split(",",$arch_hash{$b}));
		       $numb <=> $numa}  (keys(%arch_hash))){
  
  my $genelist=$arch_hash{$key};
  my $locus_number = scalar(split(",",$key));
  my $gene_number = scalar(split(",", $genelist));
  #build genelist => architecture report
  #originally due to the amount of data, we required architectures to have at least two
  #locuses, in this version we do not.
  #we now add the architecture found to the ones for this gene list.
  #the geneGroupHash{} will be used to list architectures shared by a set of genes
#  if(( ($locus)|| ($locus_number >= 2 )) && ($gene_number >= 2)){
  if($gene_number >= 2){
    $geneGroupHash{$genelist} .= "[$key]";
  }
  #prepare for printing
  $genelist =~ s/\,/\n/g;
  #  if(( ($option =~ m/locus/)|| ($locus_number >= 2 )) && ($gene_number >= 2)){
  if($gene_number >= 2){
    print "Arquitecture:$key (Locuses:$locus_number, Genes=$gene_number)\nGenes:\n$genelist\n";
    print "======================================\n";
  }
}

#now we print te report that indicates which groups of genes share which architectures
#in this version we will cross-reference the gene groups against those found in the expression
#dendograms. Only groups supported by the dendograms (that is that correspond to all genes in
#a subtree of the dendogram) are printed.

print "\n\nGENE FAMILY PATTERNS REPORT:\n===========================\n";
#the report will be sorted by the number of genes in the group
foreach my $key (sort {my $numa = scalar(split(/\,/,$a));
		    my $numb = scalar(split(/\,/,$b));
		       $numb <=> $numa}
		 (keys(%geneGroupHash))){
  #$key contains the list of genes
  #we need to sort to check agains the expression cluster hash
  my $sortedGeneGroup = join(",", sort(split(",",$key)));
  print STDERR "SORTED GENE GROUP:$sortedGeneGroup,ORIGINAL:$key\n";
  my $geneFamilyPrintFormat = "\t$sortedGeneGroup";
  $geneFamilyPrintFormat =~ s/,/\n\t/g;
  print STDERR "\t\t&&&&&&&PRINTABLE\n$geneFamilyPrintFormat\n&&&&&&&&&&&&&7\n";
  #extracth architecture list and convert to a printable format
  #in the hash each entry consists of architectures surrounded by "[" and "]".
  my $architectures = "\t".$geneGroupHash{$key};
  $architectures =~ s/\]\[/\n\t/g;
  $architectures =~ s/\[//g;
  $architectures =~ s/\]//g;
  my @archs = split("\n\t",$architectures);
  #sort alfabetically by tfbs architecture, this should naturally group some subsets
  #however we need set operations
  #$architectures = join("\n\t",(sort(@archs)));
  #old version, sorted by number of locuses, maybe we can use this order to eliminate subsets
  #we will have an array of sets.
  #first we will include the largest architecture in the array
  #next, we look at each architecture, it it is not a subset of one of the sets in the array, add it.
  $architectures = join("\n\t",(sort { scalar(split(",",$b)) <=> scalar(split(",",$a))}
  				(@archs)));
  my @setArray;
 ARCHLOOP:
  foreach my $one_architecture (sort { scalar(split(",",$b)) <=> scalar(split(",",$a))}
  				(@archs)){
    $one_architecture =~ s/\s+//g; #remove blanks
    my @tfbs_list = split(",", $one_architecture);
    my $tfbs_as_set = Set::Scalar->new();
    #
    #we are actually using multisets in a set environment, so add counters to the tfbs,this way the first and
    #the second ocurrence of a tbbs names TFBS will be addes as TFBS:0 and TFBS:1, and so forth.
    #with this trick I can successfully not report subsets.
    #for printing, the number is removed
    #
    my %tfbs_counter_hash = {};
    foreach my $one_tfbs (@tfbs_list){
      $tfbs_as_set->insert($one_tfbs.":".($tfbs_counter_hash{$one_tfbs}++));
      
    }
    #$tfbs_as_set->insert(@tfbs_list);
    print STDERR "SET OPERATIONS: checking architecture =>$tfbs_as_set<=\n";
    foreach my $set (@setArray){
      print STDERR "\tSET OPERATIONS: checking if =>$tfbs_as_set<= is contained in =>$set<=\n";
      if ($tfbs_as_set < $set){
	print STDERR "\t==>YES, go to next\n";
	next ARCHLOOP;
      }
      else{
	print STDERR "\t==>NO\n";
      }
    }
    print STDERR "SET: including set =>$tfbs_as_set<= in =>@setArray<=\n";
    push(@setArray,$tfbs_as_set);
  }
  my @genes = split(/\n\t/,$geneFamilyPrintFormat);
  foreach my $oneGene (@genes){
    $oneGene =~ s/\s+//;
    $geneFamilyHash{$oneGene} ++;
  }
  #finally, join the sets in @set_array as strings in the $architectures variable
  #this version of the program only print "gene families" that correspond to an
  #expression cluster
  if ($treefile){
    #
    #if we are considering expression data....
    #
    if ($$expression_cluster_hash{$sortedGeneGroup} == 1){
      print "-----------------------\n";
      print "GeneFamily (NEW VERSION)(",scalar(@genes)," genes):\n$geneFamilyPrintFormat\nCOMMON ARCHITECTURES(",scalar(@setArray)," architectures):\n";
      foreach my $one_arch (@setArray){
	#
	#printing the architecture, for clarity we will remove the occurence number
	#
	my $arch_string = join(",", $one_arch);
	
	$arch_string =~ s/\:\d+//g;
	print $arch_string;
      }
    }
    else{
      print STDERR "NOT COEXPRESSED GeneFamily(",scalar(@genes)," genes):\n$geneFamilyPrintFormat\nCOMMON ARCHITECTURES(",scalar(@archs)," architectures):\n$architectures\n";
    }
  }
  else {
    print "-----------------------\n";
    print "GeneFamily (NEW VERSION)(",scalar(@genes)," genes):\n$geneFamilyPrintFormat\nCOMMON ARCHITECTURES(",scalar(@setArray)," architectures):\n";
    foreach my $one_arch (@setArray){
	#
	#printing the architecture, for clarity we will remove the occurence number
	#
 	my $arch_string = join(",", $one_arch);
	
	$arch_string =~ s/\:\d+//g;
	print $arch_string;
 #     print "$one_arch\n";
    }
  }
}
print "\n\nGENE REPORT:\n=============\n";
foreach my $gene (sort {$geneFamilyHash{$b} <=> $geneFamilyHash{$a}} (keys(%geneFamilyHash))){
  print "$gene,$geneFamilyHash{$gene} families\n"
}
print "\n\nPATTERN REPORT:\n=============\n";

foreach my $pattern (sort {scalar(split("\n",$patternHash{$b})) <=> scalar(split("\n",$patternHash{$a}))}
		     (keys(%patternHash))){
  my $number = (scalar(split("\n",$patternHash{$pattern})));
  if ($number >= 2){
    print "------------------------\nPattern:$pattern\nSequences:\n$patternHash{$pattern} \n"
  }
}

print STDERR "TOTAL NUMBER OF PATTERNS: $pattern_count\n";
print STDERR "Total number of architectures:", scalar(keys(%arch_hash)), "\n";

sub subsets(){
  my $set_list = shift(@_);
  my $type = shift(@_);
  my $position = shift(@_);
  #print STDERR "ENTERING subsets($set_list,$type,$position)\n";
  #print STDERR "subsets($set_list,$type)\n";
  my $result = "";
  #we asume set is represented as a string, elements
  #separated by a colon
  #in this version locuses are defined by looking at positions, no preprocessing assumed
  my $first_tfbs_name;
  my $first_tfbs_start;
  my $first_tfbs_end;
  if ($set_list !~ m/\;/){
    #print STDERR "=>no subsets($set_list,$type)\n";
    $set_list =~ m/\s*([\w\d\W]+)\((\d+)\:(\d+)\)/; 
    my $first_tfbs_name = $1;
    my $first_tfbs_start = $2;
    my $first_tfbs_end = $3;
    if ($first_tfbs_start > $position){
      $result = $set_list;
    }
    else{
      #print STDERR "only one, overlaps, tfbs:$set_list, position=$position\n";
      $result = "";
    }
  }
  else{
    #
    #more than one TFBS
    #
    #print STDERR "more than one tfbs\n";
    my @elements = split(";", $set_list);
    my $first;
    #find the first TFBS that does not overlap previous position restriction

    my $first_tfbs_start;
    my $first_tfbs_end;
    my $first_tfbs_name;
  LOOP:
    while ($first = shift(@elements)){
      #print STDERR "cheking if $first overlaps position $position\n";
      #check fields
      $first =~ m/\s*([\w\d\W]+)\((\d+)\:(\d+)\)/;
      $first_tfbs_name = $1;
      $first_tfbs_start = $2;
      $first_tfbs_end = $3;
      #print STDERR "\t parsing $first, name=$first_tfbs_name,start=$first_tfbs_start,end=$first_tfbs_end\n";
      if ($first_tfbs_start > $position){
	#we are clear of position restrictions
	#print STDERR "\t==>no, go ahead\n";
	last LOOP;
      }
    }
    #add first tfbs to results
    $result = $first;
    if (!$locus){
      #print STDERR "NOT USING LOCUS!\n";
      my $other_elements = join(";", @elements);
      if (scalar(@elements) > 0) {
	#print STDERR "TFBSs, compute subsets first for:$other_elements, after position $first_tfbs_end\n";
	my $subset_string_with_current = subsets($other_elements,$type, $first_tfbs_end);
	my @subsets = split(";",$subset_string_with_current);
	foreach my $subset (@subsets){
	  #print STDERR "now add current TFBS ($first) to next computed subset ($subset)\n";
	  $result .= ";$first,$subset"
	}
	my $subset_string_without_current = subsets($other_elements,$type,-1);
	$result .= ";$subset_string_without_current";
	
      }
    }
    
    
    else{
      my $current_tfbs_locus = $first;
      my $last_tfbs_end = $first_tfbs_end;
    LOOP2:
      while (my $next_tfbs = shift(@elements)){
	$first =~ m/\s*[\w\d\W]+\(\d+\:\d+\)/;
	my $next_tfbs_name = $1;
	my $next_tfbs_start = $2;
	my $next_tfbs_end = $3;
	if ($next_tfbs_start > $last_tfbs_end){
	  last LOOP2;
	}
	else {
	  $current_tfbs_locus .= "/$next_tfbs";
	  $last_tfbs_end = $next_tfbs_end;
	}
      }
      #
      #took care of locus, now check if have anything left
      $result = $current_tfbs_locus;
      if (scalar(@elements) > 0) {
	my $other_elements = join(",", @elements);
	my $subset_string_with_current = subsets($other_elements,$type, $last_tfbs_end);
	my @subsets = split(";",$subset_string_with_current);
	foreach my $subset (@subsets){
	  $result .= ";$current_tfbs_locus,$subset"
	}
	my $subset_string_without_current = subsets($other_elements,$type,-1);
	$result .= ";$subset_string_without_current";
	
      }
    }
  }
  return $result;
}

sub clean(){
  my $original_subsets = shift;
  my $repeats = shift;
  my $locus = shift;
  print STDERR "entering CLEAN\n\t original string:$original_subsets)\n";
  #first remove positions, they are no longer used
  $original_subsets =~  s/\(\d+\:\d+\)//g;
  #now print
  print STDERR "\tcleaned position information:$original_subsets";
  my @architectures = split(/;/,$original_subsets);
  my %clean_set ;
  foreach my $one_arc (@architectures){
    my $final_subset;
    my @elements = split(/\,/,$one_arc);
    my @sorted_elements = sort(@elements);
    print STDERR "\tsorted information:", join(",", @sorted_elements),"\n"; 
    if ($repeats){
      #print STDERR "CONSIDERING REPEATS!!!!!Option=>$repeats<=\n"; 
      $final_subset = join(",", @sorted_elements);
    }
    else {
      #print STDERR "cleaning sorted arquitecture:",join("'", @sorted_elements),"\n";
      #in this case we remove repeats; easy since they are sorted...
      my $last_tfbs = shift(@sorted_elements);
      my @final_arc;
      push(@final_arc,$last_tfbs);
      foreach my $another_tfbs (@sorted_elements){
	print STDERR "comparing $another_tfbs with $last_tfbs\n";
	if ($another_tfbs  ne $last_tfbs){
	  print STDERR "\t OK!\n"; 
	  push(@final_arc, $another_tfbs);
	  $last_tfbs = $another_tfbs;
	  
	} 
      }
      $final_subset = join(",", @final_arc);
    }
    print STDERR "\t\tnew cleaned pattern:$final_subset\n"; 
    $clean_set{$final_subset} = 1;
  }
  my @result_list = keys(%clean_set);
  my $final_result = join(";",@result_list);
  print  STDERR "\tCLEANED SET:$final_result\n";
  return $final_result;
}

sub buildSetsFromTree{
  #to avoid using references, we will count on a global variable,
  #a hash named expressionGroups
  my %expressionGroups;
  my $tree = shift;
  chomp($tree);
  $tree =~ s/\s+//g; #remove blank spaces
  my @tree_as_array = split(undef,$tree);
  my $word = "";
  my @gene_stack;
  
   foreach my $char (@tree_as_array){
    my @aux_array;
    if ($char !~ m/[\(\)\,]/){
      #normal character, continue building word
      $word .= $char;
    }
    elsif ($char =~ /\,/){
      #comma, add word to stack
      #bit first, remove score, which are the numbers and "." after the colon
      (my $gene, my $score) = split(":", $word);
      $aux_array[0] = $gene;
      unshift( @gene_stack,@aux_array);
      $word = "";
    }
    elsif ($char =~ /\(/){
      # open parenthesis, add parethesis to stack
      $aux_array[0] = '(';
      unshift( @gene_stack,@aux_array);
    }
    elsif ($char =~  /\)/ ){
      #\t close parenthesis, build group\n
      #first add last gene name, after removing score
      (my $gene, my $score) = split(":", $word);
      $aux_array[0] = $gene;
      unshift( @gene_stack,@aux_array);
      #now process stack
      my $gene = shift(@gene_stack);
      my $gene_group = "$gene";
      while ($gene_stack[0] !~ /\(/){
	
	$gene = shift(@gene_stack);
	$gene_group .= ",$gene";
      }
      #take parenthesis out too.
      shift(@gene_stack);
      #before adding gene group to hash, sort genes by name, to avoid duplication
      my @gene_list = split(",", $gene_group);
      my $sorted_gene_group = join(",", sort(@gene_list));
      $expressionGroups{$sorted_gene_group} = 1;
      print  STDERR "EXPRESSION GROUP:$sorted_gene_group\n";
      $word = $gene_group; #next round it will add group to the stack.
    }
  }
  return \%expressionGroups;
}
