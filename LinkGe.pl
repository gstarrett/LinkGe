#!/usr/bin/perl
use warnings;
no warnings ('uninitialized', 'substr');
use strict;
use Bio::DB::Sam;
use Cwd;
use Getopt::Long;
use threads;

# Define options
my ($versionOpt, $range, $out, $help, $batch, $threaded, $fasta, $posList, $paramPath, $analysis, $chrom);
GetOptions(
			'version' => \$versionOpt,	# print program version
			'range=s' => \$range,		# use a range of positions instead of list
			'chrom=s' => \$chrom,		# chromosome or contig of interest
			'out=s' => \$out,			# user-defined outfile  prefix
			'help' => \$help,			# print usage
			'batch' => \$batch,			# process all of the bam files in a folder
			'threaded' => \$threaded,	# run each sample of a batch process on a different thread
			'FASTA' => \$fasta,			# export fasta files containing reads for each linkage polymorphism
			'list=s' => \$posList,
			'file=s' => \$paramPath,
			'Analysis=s' => \$analysis,
);

# Define version, author, contact, etc information
my $version = "0.3";
my $lastUpdate = "September 16, 2014";
my $author = "Gabriel Starrett";
my $email = "gjstarrett\@gmail.edu";
my $institution = "Originally developed at the University of Wisconsin";
my $options = <<END;
Options:

-v	--version	list current version of PEGen
-r	--range		allow an input range (e.g. 1-100)
-c	--chrom		define the chromosome or contig to analyze
-o	--out		define the prefix for the output file
-b	--batch		process all of the files in a given folder
-h	--help		show this screen
-t	--threaded	run each sample of a batch process on a different thread
-F	--FASTA		export a fasta file for each linkage pattern
-f	--file		define positions from a csv or tsv file
-l	--list		define positions from a command line list
-A	--Analysis	amino acid analysis
END
my $usage = "LinkGe\n\nVersion: $version\nLast Updated: $lastUpdate\nAuthor: $author\nEmail: $email\nInstitution: $institution\n\n$options\nUsage: perl LinkGe.pl [OPTIONS] file.bam\n\n";

# quit unless we have the correct number of command-line args
# my $num_args = $#ARGV + 1;
# if ($num_args != 2) {
#   print "Invalid number of arguments\n\n",$usage;
#   exit;
# }

my $refPath;
my @bamPath;
my $batchCount;
my @threads;
my $dir = getcwd;
my @fileNames;

#set up a subroutine for some of the options
sub options {
	if ((defined $out) && ($batch == 1)) {
		print "WARNING: cannot have user-named output file when processing in batch.\nContinuing program using default naming scheme";
	}
	if (($threaded >= 1) && ($batch != 1)) {
		print "ERROR: the threaded option only works with batch processing.\n" and die;
	}
	if ($versionOpt == 1) {
		print "PEGen version $version\n" and exit;
	}
	if ($help == 1) {
		print $usage and exit;
	}
	if ($batch == 1) {
		opendir(DIR, $ARGV[0]) or die "can't opendir $ARGV[0]: $!";
		while (defined(my $file = readdir(DIR))) {
    		if ($file =~ /bam$/) {
    			push (@bamPath, "$ARGV[0]/$file");
    		}
		}
		closedir(DIR);
	} else {
		@bamPath = $ARGV[0];
	}
}
options();
if ($threaded >= 1) {
	multithread();
} elsif (defined $analysis) {
	foreach (@bamPath) {
		my $frame = 1;
		my $bam = $_;
		$batchCount++;
#		until ($frame == 3) {
			main($bam,$frame);
			#$frame++;
#		}
	}
} else {
	foreach (@bamPath) {
	$batchCount++;
	main($_);
}
}
# The following subrouting will execute the main subroutine in a different thread for each fill being processed in batch
sub multithread {
foreach (@bamPath) {
	$batchCount++;
	push (@threads,threads->create(\&main,$_));
}
foreach (@threads) {
my @join = $_->join();
}
}

sub main {
# options();
my $startTime = time;	# start time of the program to later calculate the amount of time until completion
my @positions;
my @posHashes;
my $pairs;
my $frame = $_[1];
if ($batch == 1) {
	print "\nstarting processing for $_[0] ... $batchCount of ", scalar @bamPath,"\n";
}

# Import positions
if(defined $paramPath){
	@positions = &paramFile($paramPath);
}elsif(defined $range){
	@positions = &rangeSub($range);
}elsif(defined $posList){
	@positions = &posList($posList);
}else{
	die "Error: No positions defined\n\n$usage";
}

# import bam file
my $sam = Bio::DB::Sam->new(-bam=>"$_[0]",-fasta=> "$_[1]",-autoindex => 1) or die "Couldn't open BAM file: $!";
my @seq_ids=$sam->seq_ids;
if (defined $chrom) {
	@seq_ids=$chrom;
}
# process reference file
# print $seq_ids[0],"\n";
# print $positions[0],"\n";
# print scalar(@positions),"\n";
my $refseq;
if (defined $analysis){
	$refseq = &parseFASTA($analysis,@seq_ids,$positions[0],scalar(@positions));
}
# print "$refseq\n";
# exit;

# this runs a loop for each position listed in the param file
foreach (@positions) {
my $start = $_-1;
my $end = $_;
my ($all_trim,@all_read_names,$base_call,@read_names,@linkages,%readHash,$trim);
my $posNum = $_;

print "parsing position $_ ... " unless ($threaded >= 1);

# pull out all of the read names and base calls for reads that overlap the position and store them in an array
#my @alignments =$sam->get_features_by_location(-seq_id=>@seq_ids,-start=>$start,-end=>$end);
my $all_linkages = sub {{
                my $a = shift;
              	my @proteins;
              	my $qstart = $a->pos;
              	my $seq = $a->qseq;
               	if (defined $analysis) {
               		my $preseq = substr($seq,0,$start-$qstart);
               		@proteins = &translatorMain($preseq,$frame);
               	}
               	next if ($proteins[0] =~ /_/);
               	#print $proteins[0],"\n";
                my $base_call = substr($seq,$start-$qstart,1);
                for ($base_call) {
                	if ($base_call !~ m/A|a|T|t|C|c|G|g|U|u|N|n|G|g|P|p|V|v|L|l|I|i|M|m|F|f|Y|y|W|w|H|h|K|k|R|r|Q|q|E|e|D|d|S|s|T|t|X|x/) {
                	$base_call = "-";}
                }
                my $all_reads = $a->display_name;
            	    for ($all_reads) {
                	if (m/(.+)_./) {
                $all_trim = $1;
                $pairs++;
                	} else { $all_trim = $all_reads; }
                }
                if (not defined $range) {
            		$readHash{$all_trim} = "$base_call$posNum:";
            	} else {
            		$readHash{$all_trim} = "$base_call";
            	}
}};
$sam->pileup("$seq_ids[0]:$start-$end",$all_linkages);

print scalar (keys %readHash), " reads\n" unless ($threaded >= 1);
push @posHashes, \%readHash;
}
if ($pairs > 0) {
	print "finding read pairs that cover all positions ... " unless ($threaded >= 1);
} else { print "finding reads that cover all positions ... " unless ($threaded >= 1);
}

my $comboHash = merge(@posHashes);

my $matches = scalar keys %{$comboHash};
print $matches," matches\n" unless ($threaded >= 1);

# find read names that have the same linkage polymorphism pattern and count each occurrence and store values in a hash
print "counting reads with matching calls ... \n" unless ($threaded >= 1);
my %counts = ();
my %unique = ();
foreach my $key (sort keys %{$comboHash}) {
    my $value = ${$comboHash}{$key};
    if (not exists $counts{$value}) {
        $unique{$key} = $value;
    }
    $counts{$value}++;
};

# print frequencies to file
my $bamName;
# collect time information to provide a unique filename to each output text file
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = $hour . $minute . $second . "_" . $year . "-" . $month . "-" . $dayOfMonth;

# remove extension from bam file to name the output file
if ($_ =~ /.+\/([^\/]+)\.bam/) {
	$bamName = $1;
}
# name output file with default naming schema i.e. bamfilename_time.txt or userdefinedname.txt and open it
my $prefix = $bamName . "_". $theTime;
if (defined $out) { $prefix = $out unless ($batch == 1); } # user defined name from -o option

my %uniqueFiles;
if ($fasta == 1) {
%uniqueFiles = makeFASTA($sam,\@seq_ids,\@positions,\%{$comboHash},$prefix);
}

print "printing $prefix.txt to $dir ... \n";
open (OUT,">", "$dir/$prefix.txt");

# print header
print OUT "### LinkGe version $version\n";
if ($pairs > 1) {
	print OUT "### Total read pairs: $matches\n";
} else {
	print OUT "### Total reads: $matches\n";
}
if (defined $range){
	print OUT "### Range: $range\n"
}
if (defined $analysis){
	print OUT "### Frameshift: $frame\n"
}
my @refaa;
print OUT "id\t";
print OUT "sequence" unless (defined $analysis);
if (defined $analysis){
	@refaa = &translatorMain($refseq,$frame);
	foreach (@refaa) {
		print OUT "file name\t$_\t";
	}
}
print OUT "\treads\tfrequency\n";
foreach (keys %counts){
	print OUT "[$uniqueFiles{$_}]\t" if ($fasta == 1); # id number to easily define linkage polymorphisms also will be used to name the fasta files.
	print OUT $_,"\t" unless (defined $analysis);
	my $ntseq = $_;
	my @aaseqs;
	if (defined $analysis){
		@aaseqs = &translatorMain($ntseq,$frame);
		my $i = 0;
		foreach (@aaseqs){
			my $aaalign = &translation_align($refaa[$i],$_);
			print OUT $bamName,"\t",$aaalign,"\t";
			$i++;
		}
	}
	print OUT $counts{$_},"\t";
	printf OUT ("%.7f", $counts{$_}/(scalar keys %{$comboHash}));
	print OUT "\n";
}
close OUT;

my $endTime = time; # end time for script used to calculate total amount of run time on the following line
print "script completed in ",$endTime-$startTime,"s\n";

}

#########################
# SUBROUTINES
#########################
sub paramFile {
	my $paramPath = shift;
	my @positions;
	open (PARAM,"<", $paramPath) or die "$!";
	while (my $line = <PARAM>) {
		chomp($line);
		if ($line =~ /position_\d+(,|\t)(\d+)/g ) {
			push (@positions, $2)
		}
	}
	close PARAM;
	return @positions;
}
sub rangeSub {
	my $range = shift;
	my ($start,$end,@positions);
	#print $range."\n";
#	my $str = "ATCGCAGCTAGCTCACGACTCGCATTCTCCATACGTCACAGGTACACTACTTGACTGACGACTGATCTATACTGACTGCTGC";

	if ($range =~ /(\d+)-(\d+)/) {
		$start = $1;
		$end = $2;
	}
# 	if ($end > length($str)){
# 		die "Error: position $end is outside of sequence.";
# 	}
	my $i=$start-1;
	until ($i==$end){
		push(@positions, $i);
		$i++;
	}
# 	for my $position (@positions){
# 		print substr($str, $position, 1),"\n";
# 	}
	return @positions;
}

sub posList {
	my $arg = shift;
	my @positions = split(',',$arg);
	return @positions;
}
sub merge {
    my $first = shift;
    my @hashes = @_;
    my %result;
    KEY:
    for my $key (keys %$first) {
        my $accu = $first->{$key};
        for my $hash (@hashes) {
            next KEY unless exists $hash->{$key};
            $accu .= $hash->{$key};
        }
        $result{$key} = $accu;
    }
    return \%result;
}
sub makeFASTA {
	print "writing FASTA files to directory, $_[4] ...\n";
	my $fastaDir = "$dir/$_[4]";
	mkdir $fastaDir;
	my @positions = @{$_[2]};
	my $sam = $_[0];
	my @seq_ids = @{$_[1]};
	my %fastaHash;
	my %comboHash = %{$_[3]};
	my %unique;
	my %uniqueFiles;
	my $id = 1;
	my %fileHash;
	foreach (@positions) {
		my $start = $_-1;
		my $end = $_;
		#print Dumper(%comboHash);
		my $all_linkages = sub {
					my $a = shift;
					my $sequence = $a->qseq;
					for ($sequence) {
						if ($sequence !~ m/A|a|T|t|C|c|G|g|U|u|N|n|G|g|P|p|V|v|L|l|I|i|M|m|F|f|Y|y|W|w|H|h|K|k|R|r|Q|q|E|e|D|d|S|s|T|t|X|x/g) {
						$sequence = "-";}
					}
					my $readName = $a->display_name;
					#print ">$readName\n$sequence\n";
					$fastaHash{$readName} = "$sequence";
		};
		$sam->fetch("@seq_ids:$start-$end",$all_linkages);
	}
	foreach my $seqKey (sort keys %fastaHash) {
		my $seqValue = $fastaHash{$seqKey};
		foreach my $key (sort keys %comboHash) {
			my $value = $comboHash{$key};
			my $fileName;
			if (not exists $uniqueFiles{$value}){
				$uniqueFiles{$value} = $id;
				$fileName = $uniqueFiles{$value};
				$fileHash{$fileName} = ">$value\n";
				$id++;
			} else {$fileName = $uniqueFiles{$value}}
			if (($seqKey =~ /$key/) && (not exists $unique{$seqKey})) {
				$fileName = $uniqueFiles{$value};
				$fileHash{$fileName} .= ">$seqKey\n$seqValue\n";
				$unique{$seqKey}++;
			}
		}
	}
	#return %uniqueFiles;
	for my $fileKey (keys %fileHash) {
		open (FASTA, "> $fastaDir/$fileKey.fa");
		print FASTA $fileHash{$fileKey};
		close (FASTA);
	}
}
sub parseFASTA {
	my $fastaPath = shift;
	my $region = shift;
	my $sequence;
	my $start = shift;
	my $length = shift;
	open (REF, "<", $fastaPath) or die "Cannot open FASTA file at $fastaPath";
	while (my $line = <REF>) {
		chomp($line);
		$line =~ s/\r//g;
		if ($line !~ /^>/){
			$line =~ s/U|u/T/g;
			#print "$line\n";
			$sequence .= $line;
		}
	}
	#print $sequence,"\n";
	close (REF);
	my $subseq = substr($sequence,$start-1,$length);
	return $subseq;
}

sub translatorMain {
my $dna = shift;
$dna =~ s/-//g;
my @proteins;
	my $protein;
	for(my $i=shift; $i < (length($dna) - 2) ; $i += 3) {
		my $subseq = substr($dna,$i,3);
		#print $subseq."\n";
		$protein .= &translator($subseq);
	}
	push (@proteins, $protein);


return @proteins;
}

sub translator {
	my %daRules = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '_',    # Stop
		'TAG' => '_',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '_',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
	);

	my $codon = shift;
	if ($codon =~ /N/) {
		return "X";
	} elsif (exists $daRules{$codon}) {
		return $daRules{$codon};
	} else {die "Bad move Anton: $codon\n";}
}
sub translation_align {
	my $seq1 = shift;
	my $seq2 = shift;
	my $seqLength = length($seq1);
	#print $seqLength;
	my $finalseq;

	my $i=0;
	until ($i == $seqLength){
		my $subseq1 = substr($seq1,$i,1);
		#print "$subseq1\t";
		my $subseq2 = substr($seq2,$i,1);
		#print "$subseq2\n";
		if ($subseq2 =~ /$subseq1/) {
			$finalseq .= '.';
		} else {
			$finalseq .= "$subseq2";
		}
		$i++;
	}

	return $finalseq;
}
