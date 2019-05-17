#!/usr/bin/perl -w
#
# multiple runs: GEM algorithm
# By Jon Hallander 131126

use warnings;

die "usage: mult_run.pl <number of start populations> <name of input file> <number of runs>  <number of runs per K> <path of init file directory> <name of init file> <file extension> <save directory> <method: GEM = 1, SGEM = 2, SAGEM = 3, DAGEM = 4> <path to executable> <start number of iteration per K (to use in multiple runs, 1 for a single set of runs)> <a new starting configuration per run ( = 1) or the same configuration ( = 0)>\n" if ( $#ARGV !=11 );

chomp($K0 = $ARGV[0]);
chomp($file = $ARGV[1]);
chomp($nrun = $ARGV[2]);
chomp($nperk = $ARGV[3]);
chomp($path = $ARGV[4]);
chomp($initname = $ARGV[5]);
chomp($ext = $ARGV[6]);
chomp($dir = $ARGV[7]);
chomp($method = $ARGV[8]);
chomp($pathex = $ARGV[9]);
chomp($start = $ARGV[10]);
chomp($flag = $ARGV[11]);

open(FILEM, "<$file" ) or die "Can't open $file : $!";
chomp(@LISTM = <FILEM>);
close(FILEM);

my $Kend = $K0 + $nrun;
my $j;
my $end = $nperk + $start -1;
my $filename;

print "Enter loop...\n"; 
for ($i = $K0; $i < $Kend; $i++) {
    print "Populations, k: $i\n";
    for ($j = $start; $j <= $end; $j++){
	print "Iteration: $j\n";
	if($flag==1){
	    $filename = $path . "/" . $i . "/" . $initname . $i . "_" . $j . $ext;
	}
	else{ $filename = $path . "/" . $i . "/" . $initname . $i . "_1" . $ext; }
	print "Initialization file: $filename\n";
	my $fileout = "input_file_" . $i . "_run_" . $j . ".in";
	my $outfile = $dir . "info_run_" . $i . "_run_" . $j . ".out";
	my $outfile3 = $dir . "model_statistics_". $i . "_run_" . $j . ".out";
	my $outfile4 = $dir . "allocations_". $i . "_run_" . $j . ".out";
	open (MYFILE, '>' . $fileout) or die "Can't open $fileout : $!";

	foreach my $line (@LISTM){
	    if($line =~ m/fixed/){ print MYFILE $i . "\tFixed number of populations\n"; }
	    elsif($line =~ m/intialization/){ 
		@tmp = split('\t',$line);
		#print  $filename . "\t" . $tmp[1] . "\n"; 
		print MYFILE $filename . "\t" . $tmp[1] . "\n"; 
	    }
	    else{ print MYFILE $line . "\n"; }
	}
	close (MYFILE); 
    # add $i to input file
    # save input file
	if($method==1){
	    #my $outfile2 = $dir . "conv_alloc_". $i . ".out";
	    system("echo start GEM-algorithm...\n");
	    system("$pathex/GEM $fileout > $outfile");# run GEM
	    system("mv alloc.out $outfile4");# save output file
	    system("mv model_selection.txt $outfile3");

	}
	elsif($method==2){
	    system("echo start SGEM-algorithm...\n");
	    system("$pathex/SGEM_v1 $fileout > $outfile");# run SGEM
	   # my $outfile4 = $dir . "allocations_". $i . "_run_" . $j . ".out";
	    my $outfile6 = $dir . "mean_allocations_". $i . "_run_" . $j . ".out";
	    my $outfile7 = $dir . "mean_mixture_weights_". $i . "_run_" . $j . ".out";
	    my $outfile5 = $dir . "mixture_weights_". $i . "_run_" . $j . ".out";
	    system("mv model_selection.txt $outfile3");
	    system("mv allocations.txt $outfile4");
	    system("mv mixture_weights.txt $outfile5");
	    system("mv mean_allocation.txt $outfile6");
	    system("mv mean_mixture_proportion.txt $outfile7");
	}
	elsif($method==3){
	    system("echo start SAGEM-algorithm...\n");
	    system("$pathex/SAGEM $fileout > $outfile");# run SAGEM
	    system("mv model_selection.txt $outfile3");
	    system("mv mixture_proportion.txt $outfile4");
	}
	else{
	    system("echo start DAGEM-algorithm...\n");
	    system("$pathex/DAGEM $fileout > $outfile");# run DAGEM
	    system("mv model_selection.txt $outfile3");
	    system("mv alloc.out $outfile4");# save output file
	}
    }
	
}
