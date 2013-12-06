#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );

######################################################################################################################################################
#
#	Description
#		This is a perl script to extract line in a SAM file with the criteria specified by user.
#
#	Input
#		--samPath=						a sam file, which is coverted from a sorted bam file using sam tools, "$ samtools view -h -o out.sam in.bam"
#		--cntgRng=						the name of the contig as well as the range on the contig to be extracted, can be more than one, in the foramt of Contig:Start:End, e.g. Pf3D7_11:1256:28880; default to output all contig and all locations;
#		--strand=						strand of the read to be extracted, "+", "-", or "both"; default = "both";
#		--maxNM=						the NM attribute of sam format, number of mismatches, default = 9999, i.e. no limit; If cannot find NM, this criteria will be ignored;
#		--maxNH=						the NH attribute of sam format, number of hits, default = 9999, i.e. no limit; If cannot find NH, this criteria will be ignored;
#		--minLen=						minimum length to the read length to be extracted, default = 0, i.e. no limit;
#		--maxLen=						maximum length to the read length to be extracted, default = 9999, i.e. no limit;
#		--spliced=						"yes" or "no" to extract spliced and non-spliced reads ONLY, or "both" to extract both; default = "both";
#		--maxJunctSize=					the maximum junction size to be included, e.g. 1000, i.e. spliced read with junction size greater 1000 will be removed; default = 9999;
#		--rmRedundant=					yes or no; will remove the redundant seq if yes; default = no;
#		--rmPolyARd=					yes or no; to remove the reads with a polyA track following its last mapped genomic position; the purpose is for removing artefact polyA tailed read since the polyA tail might came from the transcript rather than polyAdenylation. The parameter was determined by downANumRegion and downANumMax option. default = no;
#		--downANumRegion=				length of the downstream (outside the read) polyA region to be checked --rmPolyARd= option; default = 10;
#		--downANumMax=					max number of As to be tolerated in the downstream (outside the read) polyA region of --downANumRegion= option default=7;
#		--upANumRegion=					length of the upstream (inside the read) polyA region to be checked --rmPolyARd= option; default = 6;
#		--upANumMax=					max number of As to be tolerated in the upstream (inside the read) polyA region of --upANumRegion= option default=4;
#		--downATrackRegion=				length of check region of downATrackMax option; default = 8;
#		--downATrackMax=				the maximum length of a consecutive track of polyA within downATrackRegion nt of downstream of the 3'End of the read, will only turn on if --rmPolyARd=yes, default = 4;
#		--refFastaPath=					path of the reference genome fasta. It is used to get the seq of the genome when rmPolyARd is turned on.
#		--countOnly=					"yes" or "no". If "yes", all other option will become irrelevant and will count the statistics only; default=no
#		--hardCodedCriteria=			"yes" or "no". If "yes", the reads will be selected by some extra hard coded criteria as coded; The aim to is give ad hoc flexibility for selection;
#		--sortSmallRNAPolyAMismatch=	"yes" or "no". If "yes", it will sort the filtered sequence into two pools, with and without polyA tail mismatch; this is designed to pick up the extra As added into Entamoeba dataset; if "yes", the default output will be split into two labelled as polyATailMsmtchPlus and polyATailMsmtchMinus, the default output wont be written; defualt = no
#		--maxHomoBasePrptn=				maximum proportion of one single base in a read, aim to prevent homopolymer reads; use 0 to switch-off; default = 0;
#		--extendEndMaxPolyA=			integer. If > 0, will extend the alignment with the genome matched As following the end of the read with maximum $extendEndMaxPolyA bases, aimed to locate the exact start of the mismatch between the polyA read and the genome seq; default = 0;
#		--dynamicDownMaxAPct			"no" or a percentage; it overrides downANumRegion and downANumMax; turn on only if rmPolyARd=yes' maxium downstream percentage of A when compare with the polyA length of the read; e.g. 60% for read with 10xA, maxium 6As in 10nt, or 60% for read with 18xA, maximum 11As in 18nt etc; will switch to 40% if polyALen of the read is < 10nt (hard-coded), more stringent;
#		--checkJunctionPolyA=			"yes" or "no", check the polyA false positves caused by junctions; if "yes", in cases which the read ended within certain bases next to the junctions, potential polyA tracks will be check on the other side of the jucntions; default = no
#		--GFFPath=						a path or "no"; if checkJunctionPolyA is yes, a GFF path must be provided; default = no
#		--outDir=						directory for output, default = "./";
#
#	Usage
#		perl SAMExtractor_v0.4.pl --samPath=/Volumes/A_MPro2TB/NGS/nov.2010/L7_noFilter/HMMSplicer/finalSAM/L7_Rahman_small.RSm.T.no.L.0.M.18.T.0.Q.18.Nr.yes_Q1_R50.edited.sorted.sam --cntgRng=DS571169:20796:28420 --cntgRng=DS571173:11508:18439 --strand=both --maxNM=3 --maxNH=5 --minLen=26 --maxLen=30 --spliced=both --countOnly=no --outDir=./test/
#
#	Assumption
#
#	History:
#		
#		V0.1
#			-debut
#
#		V0.2
#			-added --maxJunctSize= option
#
#		V0.3
#			-added --rmRedundant= option
#
#		v0.4
#			-added polyAMapping related options, --rmPolyARd=, --downANumRegion=, --downANumMax=, --refFastaPath=
#
#		v0.5
#			-added countOnly option;
#			-added hardCodedCriteria option;
#			-added sortSmallRNAPolyAMismatch option;
#
#		v0.6
#			-added maxHomoBasePrptn option
#			-added upANumRegion, upANumMax in addition to $downANumRegion, $downANumMax, aim to eliminate long polyA runs on the read lead to artifact of polyA tail;
#			-added extendEndMaxPolyA option;
#			-added samtools index and convert to bam function;
#			-added upDownPolyATrckIn10nt option;
#
#		v0.7
#			-removed upDownPolyATrckIn10nt, added upATrackRegion, upATrackMax, downATrackRegion, downATrackMax
#			-added dynamicMaxAPct option
#
#		v0.8
#			-in the remove polyRd option, added function to check false positives from junctions; the option is checkJunctionPolyA, and GFF have to be provided if checkJunctionPolyA is yes
#
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#1----------Read parameters ----------#
use vars qw ($samPath $cntgRngHsh_ref $strand $maxNM $maxNH $minLen $maxLen $spliced $maxJunctSize $rmRedundant $rmPolyARd $downANumRegion $downANumMax $upANumRegion $upANumMax $downATrackRegion $downATrackMax $refFastaPath $countOnly $hardCodedCriteria $sortSmallRNAPolyAMismatch $maxHomoBasePrptn $extendEndMaxPolyA $dynamicDownMaxAPct $checkJunctionPolyA $GFFPath $outDir $paraTag);
($samPath, $cntgRngHsh_ref, $strand, $maxNM, $maxNH, $minLen, $maxLen, $spliced, $maxJunctSize, $rmRedundant, $rmPolyARd, $downANumRegion, $downANumMax, $upANumRegion, $upANumMax, $downATrackRegion, $downATrackMax, $refFastaPath, $countOnly, $hardCodedCriteria, $sortSmallRNAPolyAMismatch, $maxHomoBasePrptn, $extendEndMaxPolyA, $dynamicDownMaxAPct, $checkJunctionPolyA, $GFFPath, $outDir, $paraTag) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

my ($rdLenHitMsmtchHsh_ref, $redundancyCountHsh_ref, $samPathHsh_ref) = readAndEditSAMOnTheFly($cntgRngHsh_ref);
runSamtoolToConvertSAMToBamAndIndex($samPathHsh_ref);
outputReadLengthHitMismatch($rdLenHitMsmtchHsh_ref);

printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$outDir = "./";
	$strand = "both";
	$maxNM = 999;
	$maxNH = 999;
	$minLen = 0;
	$maxLen = 999;
	$maxJunctSize = 999;
	$rmRedundant = "no";
	$spliced = "both";
	$rmPolyARd = "no";

	$downANumRegion = 10;
	$downANumMax = 8;
	$upANumRegion = 10;
	$upANumMax = 8;

	$downATrackRegion = 3;
	$downATrackMax = 3;

	$extendEndMaxPolyA = 0;
	
	#$refFastaPath = "/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/Pf3D7_01_v3_JAN2012_withMitoPlstd.fa";
	$refFastaPath = "/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa";
	$countOnly = "no";
	$hardCodedCriteria = "no";
	$sortSmallRNAPolyAMismatch = "no";
	$maxHomoBasePrptn = 0;
	$dynamicDownMaxAPct = "no";
	
	$checkJunctionPolyA = "no";
	$GFFPath = "no";
	
	my %cntgRngHsh;
	
	foreach my $param (@ARGV) {
		if ($param =~ m/--samPath=/) {$samPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--cntgRng=/) {
			my $cntgRng = substr ($param, index ($param, "=")+1);
			my @cntgRngSplt = split /\:/, $cntgRng;
			push @{$cntgRngHsh{$cntgRngSplt[0]}}, $cntgRngSplt[1].":".$cntgRngSplt[2];
		}
		elsif ($param =~ m/--strand=/) {$strand = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--maxNM=/) {$maxNM = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--maxNH=/) {$maxNH = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minLen=/) {$minLen = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--maxLen=/) {$maxLen = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--spliced=/) {$spliced = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--maxJunctSize=/) {$maxJunctSize = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--rmRedundant=/) {$rmRedundant = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--rmPolyARd=/) {$rmPolyARd = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--downANumRegion=/) {$downANumRegion = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--downANumMax=/) {$downANumMax = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--upANumRegion=/) {$upANumRegion = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--upANumMax=/) {$upANumMax = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--downATrackRegion=/) {$downATrackRegion = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--downATrackMax=/) {$downATrackMax = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--refFastaPath=/) {$refFastaPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--countOnly=/) {$countOnly = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--hardCodedCriteria=/) {$hardCodedCriteria = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--sortSmallRNAPolyAMismatch=/) {$sortSmallRNAPolyAMismatch = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--maxHomoBasePrptn=/) {$maxHomoBasePrptn = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--extendEndMaxPolyA=/) {$extendEndMaxPolyA = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--dynamicDownMaxAPct=/) {$dynamicDownMaxAPct = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--checkJunctionPolyA=/) {$checkJunctionPolyA = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--GFFPath=/) {$GFFPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	#---check the files
	open (TEST, "$samPath") || die "Can't open $samPath\n"; close TEST;

	if ($rmPolyARd eq "yes") {
		open (TEST, "$refFastaPath") || die "Can't open refFastaPath\n"; close TEST;
		if ($checkJunctionPolyA eq "yes") {
			open (TEST, "$GFFPath") || die "Can't open GFFPath\n"; close TEST;
		}
	}

	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	
	my $paraTag = "Len".$minLen.".".$maxLen."_NM".$maxNM."_NH".$maxNH."_rmRd".$rmRedundant;
	
	if (($downANumRegion > 20) or ($downATrackRegion > 20)) {
		die "downANumRegion and downATrackRegion must be smaller than 20\n";
	}
	if ($rmPolyARd eq "yes") {
		$paraTag = "dynPct".$dynamicDownMaxAPct."_uNum".$upANumRegion.".".$upANumMax."_dNum".$downANumRegion.".".$downANumMax."_dTrk".$downATrackRegion.".".$downATrackMax."_".$paraTag;
	}
	
	system "mkdir -p -m 777 $outDir";
	system "mkdir -p -m 777 $outDir/$paraTag";
	
	return ($samPath, \%cntgRngHsh, $strand, $maxNM, $maxNH, $minLen, $maxLen, $spliced, $maxJunctSize, $rmRedundant, $rmPolyARd, $downANumRegion, $downANumMax, $upANumRegion, $upANumMax, $downATrackRegion, $downATrackMax, $refFastaPath, $countOnly, $hardCodedCriteria, $sortSmallRNAPolyAMismatch, $maxHomoBasePrptn, $extendEndMaxPolyA, $dynamicDownMaxAPct, $checkJunctionPolyA, $GFFPath, $outDir, $paraTag);
}
########################################################################## readAndEditSAMOnTheFly
sub readAndEditSAMOnTheFly {

	my %cntgRngHsh = %{$_[0]};

	#---read the SAM
	open (INSAM, $samPath);
	open (TMPLOG, ">$outDir/tmp.log.txt");

	#---define the variables
	my $lastCntg = "intiation";
	my @samPathSplt = split /\//, $samPath;
	my %samPathHsh;
	if ($countOnly ne "yes") {
		if ($sortSmallRNAPolyAMismatch eq "no") {
			my $outSAMPath = "$outDir/$paraTag/extracted.$samPathSplt[-1]";
			$samPathHsh{"outSAMPath"} = $outSAMPath;
			open OUTSAM, ">$outSAMPath";
			if ($extendEndMaxPolyA > 0) {
				my $extendEndMaxPolyAPath = "$outDir/$paraTag/extendEndMaxPolyA.all.$samPathSplt[-1]";
				$samPathHsh{"extendEndMaxPolyAPath"} = $extendEndMaxPolyAPath;
				open EXTENDPOLYA, ">$extendEndMaxPolyAPath";
			}
		} elsif ($sortSmallRNAPolyAMismatch eq "yes") {
			my $polyATailMsmtchPlusPath = "$outDir/$paraTag/polyATailMsmtchPlus.all.$samPathSplt[-1]";
			my $polyATailMsmtchMinusPath = "$outDir/$paraTag/polyATailMsmtchMinus.$samPathSplt[-1]";
			$samPathHsh{"polyATailMsmtchPlusPath"} = $polyATailMsmtchPlusPath;
			$samPathHsh{"polyATailMsmtchMinusPath"} = $polyATailMsmtchMinusPath;
			open ATAILPLUSSAM, ">$polyATailMsmtchPlusPath";
			open ATAILMINUSSAM, ">$polyATailMsmtchMinusPath";
		} else {
			die "invalid sortSmallRNAPolyAMismatch option\n";
		}
	}
	
	#---get the total line number and determine the interval size
	my ($fileTotalLineNum,  $intervalSize)= checkFileSizeAndDefineIntervalSize($samPath, 100000);

	#---define the start time and counters
	my $intervalStart = time();
	my $lineProc = my $progCount = my $totalTimeSpent = 0;

	#---define the SAM Table
	my $SAMFlagTableHsh_ref = defineSAMFlagTable();
	my %SAMFlagTableHsh = %{$SAMFlagTableHsh_ref};

	my $extractCtngNum = keys %cntgRngHsh;
	my $extractAllCntg = "no";
	$extractAllCntg = "yes" if ($extractCtngNum == 0);

	my %cntgPosLenHsh;
	
	#---get the refFastaSeq if it rm polyARd is on
	my $refFastaHsh_ref;
	my %refFastaHsh;
	my $twoWaysJunctIndxHsh_ref;
	my %twoWaysJunctIndxHsh;
	
 	if ($rmPolyARd eq "yes") {
		$refFastaHsh_ref = readMultiFasta($refFastaPath);
		%refFastaHsh = %{$refFastaHsh_ref};
		if ($checkJunctionPolyA eq "yes") {
			my ($nameByGeneHsh_ref, $strndByGeneHsh_ref, $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $exonNumByCntgHsh_ref, $geneExonLenHsh_ref, $geneCtgryHsh_ref, $ctgryReadCountHsh_ref, $intronRngByGeneHsh_ref) = readGff($GFFPath);
			$twoWaysJunctIndxHsh_ref = generateTwoWaysJunctionIndex($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref);
			%twoWaysJunctIndxHsh = %{$twoWaysJunctIndxHsh_ref};
		}
	}
	
	my %rdLenHitMsmtchHsh; #--- to store the info of length, mismatch and number of hit count;
	my %redundancyCountHsh; #--- to store the info redundancy;

	my (%rdLenPolyATailCountHsh, %rdLenNoPolyATailCountHsh, %polyAPreCalSeqHsh);
	
	print "Start process the alignments.\n";
	printProgressScale("Extracting the SAM file", 50);

	#---Start reading samPath
	while (my $theLine = <INSAM>) {

		# HTt83343643	0	EhR1rDNA	11053	255	22M611N16M	*	0	0	ATTCCCACTGTCCCTATCTGCATTTCAAGCAGAATTGA	7:767)8777D=DDDDADCD=@6@?A?AAA@B@@?A>A	XA:i:0	XS:A:+
		
		#---report the progress
		chomp $theLine;
		
		#---check header, if yes, print it out and next line
		if ($theLine =~ m/^\@/) {#---skip the commment and info lines
			if ($countOnly ne "yes") {
				if ($sortSmallRNAPolyAMismatch eq "no") {
					print OUTSAM $theLine."\n";
				} elsif ($sortSmallRNAPolyAMismatch eq "yes") {
					print ATAILPLUSSAM $theLine."\n";
					print ATAILMINUSSAM $theLine."\n";
				}
			}
			next;
		}
		
		#---check progress
		$lineProc++; $progCount++;
		if ($progCount >= $intervalSize) {
			($progCount, $intervalStart) = reportProgress($progCount, $lineProc, $intervalSize, $fileTotalLineNum, $intervalStart);
		}
		
		#---split the line if not header
		my @theLineSplt = split /\t/, $theLine;

		#---get the length, NH and NM
		my $readStart = $theLineSplt[3];
		my $cigarStr = $theLineSplt[5];
		my $alignSeq = $theLineSplt[9];
		my $alignQual = $theLineSplt[10];
		my $readLength = length $alignSeq;

		#---check macthed position
		my $genomicLength = 0;
		my @matchLengthAry;
		while ($cigarStr =~ /(\d+)M/g) {
			$genomicLength += $1;
			push @matchLengthAry, $1;
		}

		#---check splicing
		my $readSpliced = "no";
		$readSpliced = "yes" if ($cigarStr =~ m/N/);

		#---check splicing junction size
		my $longestJunction = 0;
		while ($cigarStr =~ /(\d+)N/g) {
			$genomicLength += $1;
			$longestJunction = $1 if $1 > $longestJunction;
		}
		
		my $NH = my $NM = my $MD = -1;
		foreach my $attribute (@theLineSplt) {
			if ($attribute =~ m/NH\:i\:/) {
				$NH = $attribute;
				$NH =~ s/NH\:i\://;
			}
			if ($attribute =~ m/NM\:i\:/) {
				$NM = $attribute;
				$NM =~ s/NM\:i\://;
			}
			if ($attribute =~ m/MD\:Z\:/) {
				$MD = $attribute;
				$MD =~ s/MD\:Z\://;
			}
		}
		
		#--skip if not on cntg rng;
		my $curntCntg = $theLineSplt[2];
		
		next if $curntCntg eq "*"; #---unaligned

		#--skip if out of cntg rng;
		my $inRng = "no";
		foreach my $range (@{$cntgRngHsh{$curntCntg}}) {
			my @rangeSplt = split /\:/, $range;
			my $rngStart = $rangeSplt[0];
			my $rngEnd = $rangeSplt[1];
			$inRng = "yes" if (($readStart >= $rngStart) and ($readStart <= $rngEnd));
		}
		
		#---store all SAM bits in a Hsh
		my $SAMFlag = $theLineSplt[1];
		my $SAMBitStr = $SAMFlagTableHsh{$SAMFlag};
		my @SAMBitAry = split /\+/, $SAMBitStr;
		my %SAMBitHsh;
		foreach my $SAMBit (@SAMBitAry) {$SAMBitHsh{$SAMBit}++;}
		
		#---check strand
		my $readStrand;
		if (not exists $SAMBitHsh{16}) {$readStrand = "+";} 
		else {$readStrand = "-";}

		#--skip if redundant
		my $redundant = "no";
		if (exists ${${${$cntgPosLenHsh{$curntCntg}}{$readStart}}{$readStrand}}{$readLength}) {
			$redundant = "yes";
		}
		
		#---skip if polyA is at 5'end
		my $polyA = "no";
		
		if ($rmPolyARd eq "yes") {
			
			#---check whether sequence already is present
			my $read3End = $readStart + $genomicLength - 1;
			$read3End = $readStart if ($readStrand eq "-");
			
			if (not exists ${${$polyAPreCalSeqHsh{$curntCntg}}{$read3End}}{$readStrand}) {#---if the sequence are stored
				
				my ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean);
				
				#------define the start and end region then get the sequences
				my $cntgSeq = $refFastaHsh{$curntCntg};
				my $down20ntStart = my $up20ntStart = 0;
				my $down20ntGenSeq = my $up20ntGenSeq = "";
				if ($readStrand eq "+") {
					$down20ntStart = $read3End;
					$up20ntStart = $read3End - 20 - 1;
					$down20ntGenSeq = substr $cntgSeq, $down20ntStart, 20;
					$up20ntGenSeq = substr $cntgSeq, $up20ntStart, 20;
				} elsif ($readStrand eq "-") {
					$down20ntStart = $readStart - 20 - 1;
					$up20ntStart = $readStart;
					$down20ntGenSeq = substr $cntgSeq, $down20ntStart, 20;
					$up20ntGenSeq = substr $cntgSeq, $up20ntStart, 20;
					$down20ntGenSeq = reverse($down20ntGenSeq); $down20ntGenSeq =~ tr/ACGTacgt/TGCAtgca/;
					$up20ntGenSeq = reverse($up20ntGenSeq); $up20ntGenSeq =~ tr/ACGTacgt/TGCAtgca/;
				} else {
					die "undefined read strand\n";
				}

				my $downRegNumCheckSeq = substr $down20ntGenSeq, 0, $downANumRegion;
				my $downRegTrackCheckSeq = substr $down20ntGenSeq, 0, $downATrackRegion;
				my $upRegNumCheckSeq = substr $up20ntGenSeq, 0, $upANumRegion;
				
				#----check whether the read ends close to junction
				my $nextToJunct = 'no';
				my $down20ntcDNA = '';
				
				if ($checkJunctionPolyA eq "yes") {
					my @searchRngAry = ($read3End..($read3End+20));
					@searchRngAry = (($read3End-20)..$read3End) if ($readStrand eq "-");
					my $junctStart = my $junctEnd = 0;
					foreach my $srchPos (@searchRngAry) {#---check if the junction is around
						if (exists ${${$twoWaysJunctIndxHsh{$curntCntg}}{$readStrand}}{$srchPos}) {
							$junctStart = $srchPos;
							$junctEnd = ${${$twoWaysJunctIndxHsh{$curntCntg}}{$readStrand}}{$srchPos};
							$nextToJunct = "yes";
							last;
						}
					}

					#---get the "spliced" sequences, head half and the tail half
					if ($nextToJunct eq "yes") {
						if ($readStrand eq "+") {
							my $leftHalfLength = $junctStart - $read3End - 1;
							my $leftHalfSeq = substr $cntgSeq, $read3End, $leftHalfLength;
							my $rightHalfLength = 20 - $leftHalfLength;
							my $rightHalfSeq = substr $cntgSeq, $junctEnd, $rightHalfLength;
							$down20ntcDNA = $leftHalfSeq.$rightHalfSeq;
			
						} elsif ($readStrand eq "-") {
							my $rightHalfLength = $read3End - $junctStart - 1;
							my $rightHalfSeq = substr $cntgSeq, $junctStart, $rightHalfLength;
							my $leftHalfLength = 20 - $rightHalfLength;
							my $leftHalfSeq = substr $cntgSeq, ($junctEnd - $leftHalfLength-1) , $leftHalfLength;
							$down20ntcDNA = $leftHalfSeq.$rightHalfSeq;
							$down20ntcDNA = reverse($down20ntcDNA); $down20ntcDNA =~ tr/ACGTacgt/TGCAtgca/;
							
						} else {
							die "undefined read strand\n";
						}
					}
				}

				$withinJunct20ntBoolean = $nextToJunct;
				$genomicCheckSeq20nt = $down20ntGenSeq;
				$cDNACheckSeq20nt = $down20ntcDNA;
				$downATrackBoolean = 'no';
				$downANumBoolean = 'no';
				$upANumBoolean = 'no';
				$downANumBoolean = 'yes' if ((($downRegNumCheckSeq =~ tr/A//) > $downANumMax) and ($downANumMax > 0));
				$upANumBoolean = 'yes' if ((($upRegNumCheckSeq =~ tr/A//) > $upANumMax) and ($upANumMax > 0));
				$downATrackBoolean = 'yes' if (($downRegTrackCheckSeq =~ m/([A]{$downATrackMax,})/) and ($downATrackMax > 0));
				
				@{${${$polyAPreCalSeqHsh{$curntCntg}}{$read3End}}{$readStrand}} = ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean);
				print TMPLOG join "", ((join "\t", (($curntCntg.":".$read3End.":".$readStrand), @{${${$polyAPreCalSeqHsh{$curntCntg}}{$read3End}}{$readStrand}})), "\n");
			}
			
			#----get the polyA Seq data
			my ($withinJunct20ntBoolean, $genomicCheckSeq20nt, $cDNACheckSeq20nt, $downATrackBoolean, $downANumBoolean, $upANumBoolean) = @{${${$polyAPreCalSeqHsh{$curntCntg}}{$read3End}}{$readStrand}};
			
			#----check the fixed genomic boolean first
			if (($downATrackBoolean eq "yes") or ($downANumBoolean eq "yes") or ($upANumBoolean eq "yes")) {
				$polyA = "yes";
			
			} else {#----check the dynamic region
				
				#------get the polyA tail Length and define the downADynRegion and downADynMax;
				my $downADynRegion = $downANumRegion;
				my $downADynMax = $downANumMax;
				my $rdName = $theLineSplt[0];
				if ($rdName =~ m/\[Ax(\d+)\]/) {#---defined polyA
					my $rdPolyANum = $1;
					$rdPolyANum = 20 if ($rdPolyANum > 20);
					$downADynRegion = $rdPolyANum;
					if ($rdPolyANum > 10) {
						$downADynMax = sprintf "%.0f", $downADynRegion*($dynamicDownMaxAPct/100);
					} elsif (($rdPolyANum > 5) and ($rdPolyANum <=10)) {
						$downADynMax = sprintf "%.0f", $downADynRegion*(40/100);
					} else {#---<=5
						$downADynMax = sprintf "%.0f", $downADynRegion*(20/100) ;
					}
				} else { #---undefined polyA
					$polyA = "yes";
				}
				
				#----check the dyn seq region
				my $downRegDynCheckSeq = substr $genomicCheckSeq20nt, 0, $downADynRegion;
				$polyA = "yes" if ((($downRegDynCheckSeq =~ tr/A//) > $downADynMax) and ($downADynMax > 0));
				
				#----check the cDNA is close to junction
				if ($withinJunct20ntBoolean eq "yes") {
					my $downRegNumCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downANumRegion;
					my $downRegTrackCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downATrackRegion;
					my $downRegDynCheckcDNASeq = substr $cDNACheckSeq20nt, 0, $downADynRegion;
					if (((($downRegDynCheckcDNASeq =~ tr/A//) > $downADynMax) and ($downADynMax > 0))
					or ((($downRegNumCheckcDNASeq =~ tr/A//) > $downANumMax) and ($downANumMax > 0))
					or (($downRegTrackCheckcDNASeq =~ m/([A]{$downATrackMax,})/) and ($downATrackMax > 0))) {
						$polyA = "yes";
					}
				}
			}
		}
		
		#---record the read info
		${$redundancyCountHsh{$readLength}}{"total"}++;
		${$redundancyCountHsh{$readLength}}{"redundant"}++ if ($redundant eq "yes");

		#--- store the info
		${${$rdLenHitMsmtchHsh{$readLength}}{$NM}}{$NH}++;

		#---check maxHomoBasePrptn
		my $exceedMaxHomoBasePrptn = "no";
		if ($maxHomoBasePrptn != 0) {
			my $maxCount = int ($readLength*$maxHomoBasePrptn);
			my $ACount = ($alignSeq =~ tr/Aa//);
			my $TCount = ($alignSeq =~ tr/Tt//);
			my $GCount = ($alignSeq =~ tr/Gg//);
			my $CCount = ($alignSeq =~ tr/Cc//);
			$exceedMaxHomoBasePrptn = "yes" if (($ACount > $maxCount) or ($TCount > $maxCount) or ($GCount > $maxCount) or ($CCount > $maxCount));
		}

		next if ((not exists $cntgRngHsh{$curntCntg}) and ($extractAllCntg eq "no"));
		next if (($inRng eq "no") and ($extractAllCntg eq "no"));
		next if (($readLength < $minLen) or ($readLength > $maxLen));
		next if (($NM > $maxNM) or ($NH > $maxNH));
		next if (($readSpliced ne $spliced) and ($spliced ne "both"));
		next if ($longestJunction >= $maxJunctSize);
		next if (($readStrand ne $strand) and ($strand ne "both"));
		next if (($rmRedundant eq "yes") and ($redundant eq "yes"));
		next if (($rmPolyARd eq "yes") and ($polyA eq "yes"));
		next if (exists $SAMBitHsh{4});
		next if ($exceedMaxHomoBasePrptn eq "yes");
		
		if ($hardCodedCriteria eq "yes") {
			#----check end mismatch
			my @MDAry = split /\D/, $MD;
			@MDAry = reverse @MDAry if ($readStrand eq "+");
			my $relativePos = 0;
			my $mismatchWithinLast5nt = "no";
			foreach my $matchLength (@MDAry) {
				next if length $matchLength < 1;
				$relativePos += $matchLength;
				$mismatchWithinLast5nt = "yes" if $relativePos < 5;
			}
			
			next if ($mismatchWithinLast5nt eq "yes");
			
			next if (($matchLengthAry[-1] < 5) or ($matchLengthAry[0] < 5)); #----end match block < 5
			
			#------To remove short read mismatch 
			next if (($readLength <= 20) and ($NM > 0)); #---all shorter than 20nt reads have to be 100% match
			next if (($readLength <= 21) and ($NM > 1)); #---all shorter than 21nt reads can only have 1 mismatch
			next if (($readLength <= 22) and ($NM > 2)); #---all shorter than 22nt reads can only have 2 mismatch
			next if (($NH > 1) and ($NM > 0)); #---multiple hit doesnt allow to have mismatch
			#next if (($readLength < 20) and ($NH > 1)); #---all shorter than 20nt reads have to be unique
			#next if (($readLength < 30) and ($NH > 10)); #---all shorter than 30nt reads less than 10 hits
		}
	
		if ($sortSmallRNAPolyAMismatch ne "yes") {
			
			if ($countOnly ne "yes") {
				
				if ($extendEndMaxPolyA > 0) {
					my $extended = "no";
					if ($readStrand eq "+") {
						my $downRegNumCheckStart = $readStart + $readLength - 1;
						my $downRegNumCheckSeq = substr $refFastaHsh{$curntCntg}, $downRegNumCheckStart, $extendEndMaxPolyA;
						my $extendANum = 0;
						my @downRegNumCheckSeqAry = split //, $downRegNumCheckSeq;
						foreach my $base (@downRegNumCheckSeqAry) {
							last if (($base ne "A") and ($base ne "a"));
							$extendANum++;
						}
						
						if ($extendANum > 0) {
							$extended  = "yes";
							my $extendedSeq = "";
							my $extendedQual = "";
							foreach (1..$extendANum) {
								$extendedSeq .= "A";
								$extendedQual .= "I";
							}
							
							$extendedQual = $alignQual.$extendedQual;
							$extendedSeq = $alignSeq.$extendedSeq;
							my $extendedLength = $readLength + $extendANum;
							my $extendedCigarStr = $extendedLength."M";
							my $extendedReadStart = $readStart;
							$theLineSplt[3] = $extendedReadStart;
							$theLineSplt[5] = $extendedCigarStr;
							$theLineSplt[9] = $extendedSeq;
							$theLineSplt[10] = $extendedQual;
						}
						
					} elsif ($readStrand eq "-") {	
					
						my $downRegNumCheckStart = $readStart - $extendEndMaxPolyA - 1;
						my $downRegNumCheckSeq = substr $refFastaHsh{$curntCntg}, $downRegNumCheckStart, $extendEndMaxPolyA;
						my $extendANum = 0;
						my @downRegNumCheckSeqAry = split //, $downRegNumCheckSeq;
						
						@downRegNumCheckSeqAry = reverse(@downRegNumCheckSeqAry);
						
						#---search backwards
						foreach my $base (@downRegNumCheckSeqAry) {
							last if (($base ne "T") and ($base ne "t"));
							$extendANum++;
						}
						
						if ($extendANum > 0) {
							$extended  = "yes";
							my $extendedSeq = "";
							my $extendedQual = "";
							foreach (1..$extendANum) {
								$extendedSeq .= "T";
								$extendedQual .= "I";
							}
							
							$extendedQual = $extendedQual.$alignQual;
							$extendedSeq = $extendedSeq.$alignSeq;
							my $extendedLength = $readLength + $extendANum;
							my $extendedCigarStr = $extendedLength."M";
							my $extendedReadStart = $readStart - $extendANum;
							$theLineSplt[3] = $extendedReadStart;
							$theLineSplt[5] = $extendedCigarStr;
							$theLineSplt[9] = $extendedSeq;
							$theLineSplt[10] = $extendedQual;
						}
					}
					
					my $extendedLine = join "\t", @theLineSplt;#----the contents are extended, if any
					print EXTENDPOLYA $extendedLine."\n";

				}
			
				print OUTSAM $theLine."\n";
			}
			
		} elsif ($sortSmallRNAPolyAMismatch eq "yes") {
		
			#------To remove any read that is not having A mismatches at the end, aimed to pickup the polyA small RNAs
			my $validPolyATailMismatch = "no";

			if ($NM > 1) {
			
				#---valid if has a minATailLen lenght of polyA tail with minEndAMismatchPrptn proportion of it is mismatched, e.g. polyA tail of 2 with 1 mistach,  polyA tail of 1 with 1 mistach, polyA tail of 5 with 1 mistach
				my $minEndAMismatchPrptn = 0.2;
				my $minATailLen = 1;
				
				my $rdSeq = $alignSeq;
				if ($readStrand eq "-") {
					$rdSeq = reverse $rdSeq;
					$rdSeq =~ tr/ACGTacgt/TGCAtgca/;
				}
				
				my $ATailLen = 0;
				$ATailLen = length $1 if ($rdSeq =~ m/([A]{1,})$/);
				
				if ($ATailLen >= $minATailLen) {
					#actual read 	AGAAATTAGTAGAAAAGAAAAAAAA   MD:Z:0T18G4T0   NM:i:3  NH:i:1
					#genome seq 	TGAAATTAGTAGAAAAGAAGAAAAT	
					my @MDSeqSplt = split /A|T|G|C|N/, $MD; #---0, 18, 4, 0;
					my @MDNumSplt = split /\d+/, $MD; #---T, G, T;
					my %mismtchPosHsh;
					my $curntPos = 0;
					my $mismatchNum = @MDNumSplt;
					my @alignSeqSplt = split //, $alignSeq; 
					for my $index (0..($mismatchNum-1)) {
						my $numResidue = $MDSeqSplt[$index];
						my $genomeResSeq = $MDNumSplt[$index];
						$curntPos += $numResidue;
						my $alignResSeq = $alignSeqSplt[$curntPos];
						@{$mismtchPosHsh{$curntPos}} = ($genomeResSeq, $alignResSeq); #--- 0->T, 19->G, 24->T
						$curntPos++;
					}
					my @searchPosAry = ();
					if ($readStrand eq "+") {
						@searchPosAry = sort {$b <=> $a} (($readLength - 1 - $ATailLen - 1) .. ($readLength - 1)); #---search from end
					} elsif ($readStrand eq "-") {
						@searchPosAry = sort {$a <=> $b} (0..($ATailLen - 1)); #---search from beginning
					} else {
						die "Strand not specified.\n";
					}

					my $endAMismatchNum = 0;
					for my $pos (@searchPosAry) {
						$endAMismatchNum++ if (exists $mismtchPosHsh{$pos});
					}
							
					my $endAMismatchPrptn = $endAMismatchNum/$ATailLen;

					if ($endAMismatchPrptn >= $minEndAMismatchPrptn) {
						${${$rdLenPolyATailCountHsh{$readLength}}{$ATailLen}}{$endAMismatchNum}++;
						$validPolyATailMismatch = "yes";
					} 
				}
			}

			if ($validPolyATailMismatch eq "no") {
				$rdLenNoPolyATailCountHsh{$readLength}++;
			}
			
			if ($countOnly ne "yes") {
				if ($validPolyATailMismatch eq "yes") {
					print ATAILPLUSSAM $theLine."\n";
				} elsif ($countOnly ne "yes") {
					print ATAILMINUSSAM $theLine."\n";
				}
			}
		}
	
		#---store redudancy
		${${${$cntgPosLenHsh{$curntCntg}}{$readStart}}{$readStrand}}{$readLength}++ if ($rmRedundant eq "yes");
	}

	print  "\n";
	
	printMismatchPolyATailInfo(\%rdLenPolyATailCountHsh, \%rdLenNoPolyATailCountHsh) if ($sortSmallRNAPolyAMismatch eq "yes");
	
	close INSAM;
	
	if ($countOnly ne "yes") {
		if ($sortSmallRNAPolyAMismatch eq "no") {
			close OUTSAM;
			if ($extendEndMaxPolyA > 0) {
				close EXTENDPOLYA;
			}
		} elsif ($sortSmallRNAPolyAMismatch eq "yes") {
			close ATAILPLUSSAM;
			close ATAILMINUSSAM;
		}
	}
	
	return \%rdLenHitMsmtchHsh, \%redundancyCountHsh, \%samPathHsh;
}
########################################################################## defineSAMFlagTable
sub defineSAMFlagTable {
	
#
#copied from http://bioinformatics.bc.edu/chuanglab/wiki/index.php/SAM_pairwise_flag_translation_table
#
#0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair 
#0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)  
#0x0004 4 the query sequence itself is unmapped 
#0x0008 8 the mate is unmapped  
#0x0010 16 strand of the query (0 for forward; 1 for reverse strand) 
#0x0020 32 strand of the mate  
#0x0040 64 the read is the ﬁrst read in a pair  
#0x0080 128 the read is the second read in a pair 
#0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records) 
	
	my %SAMFlagTableHsh = (

		0 => "0",
		1 => "1",
		2 => "2",
		3 => "1+2",
		4 => "0+4",
		5 => "1+4",
		6 => "0+2+4",
		7 => "1+2+4",
		8 => "0+8",
		9 => "1+8",
		10 => "0+2+8",
		11 => "1+2+8",
		12 => "0+4+8",
		13 => "1+4+8",
		14 => "0+2+4+8",
		15 => "1+2+4+8",
		16 => "0+16",
		17 => "1+16",
		18 => "0+2+16",
		19 => "1+2+16",
		20 => "0+4+16",
		21 => "1+4+16",
		22 => "0+2+4+16",
		23 => "1+2+4+16",
		24 => "0+8+16",
		25 => "1+8+16",
		26 => "0+2+8+16",
		27 => "1+2+8+16",
		28 => "0+4+8+16",
		29 => "1+4+8+16",
		30 => "0+2+4+8+16",
		31 => "1+2+4+8+16",
		32 => "0+32",
		33 => "1+32",
		34 => "0+2+32",
		35 => "1+2+32",
		36 => "0+4+32",
		37 => "1+4+32",
		38 => "0+2+4+32",
		39 => "1+2+4+32",
		40 => "0+8+32",
		41 => "1+8+32",
		42 => "0+2+8+32",
		43 => "1+2+8+32",
		44 => "0+4+8+32",
		45 => "1+4+8+32",
		46 => "0+2+4+8+32",
		47 => "1+2+4+8+32",
		48 => "0+16+32",
		49 => "1+16+32",
		50 => "0+2+16+32",
		51 => "1+2+16+32",
		52 => "0+4+16+32",
		53 => "1+4+16+32",
		54 => "0+2+4+16+32",
		55 => "1+2+4+16+32",
		56 => "0+8+16+32",
		57 => "1+8+16+32",
		58 => "0+2+8+16+32",
		59 => "1+2+8+16+32",
		60 => "0+4+8+16+32",
		61 => "1+4+8+16+32",
		62 => "0+2+4+8+16+32",
		63 => "1+2+4+8+16+32",
		64 => "0+64",
		65 => "1+64",
		66 => "0+2+64",
		67 => "1+2+64",
		68 => "0+4+64",
		69 => "1+4+64",
		70 => "0+2+4+64",
		71 => "1+2+4+64",
		72 => "0+8+64",
		73 => "1+8+64",
		74 => "0+2+8+64",
		75 => "1+2+8+64",
		76 => "0+4+8+64",
		77 => "1+4+8+64",
		78 => "0+2+4+8+64",
		79 => "1+2+4+8+64",
		80 => "0+16+64",
		81 => "1+16+64",
		82 => "0+2+16+64",
		83 => "1+2+16+64",
		84 => "0+4+16+64",
		85 => "1+4+16+64",
		86 => "0+2+4+16+64",
		87 => "1+2+4+16+64",
		88 => "0+8+16+64",
		89 => "1+8+16+64",
		90 => "0+2+8+16+64",
		91 => "1+2+8+16+64",
		92 => "0+4+8+16+64",
		93 => "1+4+8+16+64",
		94 => "0+2+4+8+16+64",
		95 => "1+2+4+8+16+64",
		96 => "0+32+64",
		97 => "1+32+64",
		98 => "0+2+32+64",
		99 => "1+2+32+64",
		100 => "0+4+32+64",
		101 => "1+4+32+64",
		102 => "0+2+4+32+64",
		103 => "1+2+4+32+64",
		104 => "0+8+32+64",
		105 => "1+8+32+64",
		106 => "0+2+8+32+64",
		107 => "1+2+8+32+64",
		108 => "0+4+8+32+64",
		109 => "1+4+8+32+64",
		110 => "0+2+4+8+32+64",
		111 => "1+2+4+8+32+64",
		112 => "0+16+32+64",
		113 => "1+16+32+64",
		114 => "0+2+16+32+64",
		115 => "1+2+16+32+64",
		116 => "0+4+16+32+64",
		117 => "1+4+16+32+64",
		118 => "0+2+4+16+32+64",
		119 => "1+2+4+16+32+64",
		120 => "0+8+16+32+64",
		121 => "1+8+16+32+64",
		122 => "0+2+8+16+32+64",
		123 => "1+2+8+16+32+64",
		124 => "0+4+8+16+32+64",
		125 => "1+4+8+16+32+64",
		126 => "0+2+4+8+16+32+64",
		127 => "1+2+4+8+16+32+64",
		128 => "0+128",
		129 => "1+128",
		130 => "0+2+128",
		131 => "1+2+128",
		132 => "0+4+128",
		133 => "1+4+128",
		134 => "0+2+4+128",
		135 => "1+2+4+128",
		136 => "0+8+128",
		137 => "1+8+128",
		138 => "0+2+8+128",
		139 => "1+2+8+128",
		140 => "0+4+8+128",
		141 => "1+4+8+128",
		142 => "0+2+4+8+128",
		143 => "1+2+4+8+128",
		144 => "0+16+128",
		145 => "1+16+128",
		146 => "0+2+16+128",
		147 => "1+2+16+128",
		148 => "0+4+16+128",
		149 => "1+4+16+128",
		150 => "0+2+4+16+128",
		151 => "1+2+4+16+128",
		152 => "0+8+16+128",
		153 => "1+8+16+128",
		154 => "0+2+8+16+128",
		155 => "1+2+8+16+128",
		156 => "0+4+8+16+128",
		157 => "1+4+8+16+128",
		158 => "0+2+4+8+16+128",
		159 => "1+2+4+8+16+128",
		160 => "0+32+128",
		161 => "1+32+128",
		162 => "0+2+32+128",
		163 => "1+2+32+128",
		164 => "0+4+32+128",
		165 => "1+4+32+128",
		166 => "0+2+4+32+128",
		167 => "1+2+4+32+128",
		168 => "0+8+32+128",
		169 => "1+8+32+128",
		170 => "0+2+8+32+128",
		171 => "1+2+8+32+128",
		172 => "0+4+8+32+128",
		173 => "1+4+8+32+128",
		174 => "0+2+4+8+32+128",
		175 => "1+2+4+8+32+128",
		176 => "0+16+32+128",
		177 => "1+16+32+128",
		178 => "0+2+16+32+128",
		179 => "1+2+16+32+128",
		180 => "0+4+16+32+128",
		181 => "1+4+16+32+128",
		182 => "0+2+4+16+32+128",
		183 => "1+2+4+16+32+128",
		184 => "0+8+16+32+128",
		185 => "1+8+16+32+128",
		186 => "0+2+8+16+32+128",
		187 => "1+2+8+16+32+128",
		188 => "0+4+8+16+32+128",
		189 => "1+4+8+16+32+128",
		190 => "0+2+4+8+16+32+128",
		191 => "1+2+4+8+16+32+128",
		192 => "0+64+128",
		193 => "1+64+128",
		194 => "0+2+64+128",
		195 => "1+2+64+128",
		196 => "0+4+64+128",
		197 => "1+4+64+128",
		198 => "0+2+4+64+128",
		199 => "1+2+4+64+128",
		200 => "0+8+64+128",
		201 => "1+8+64+128",
		202 => "0+2+8+64+128",
		203 => "1+2+8+64+128",
		204 => "0+4+8+64+128",
		205 => "1+4+8+64+128",
		206 => "0+2+4+8+64+128",
		207 => "1+2+4+8+64+128",
		208 => "0+16+64+128",
		209 => "1+16+64+128",
		210 => "0+2+16+64+128",
		211 => "1+2+16+64+128",
		212 => "0+4+16+64+128",
		213 => "1+4+16+64+128",
		214 => "0+2+4+16+64+128",
		215 => "1+2+4+16+64+128",
		216 => "0+8+16+64+128",
		217 => "1+8+16+64+128",
		218 => "0+2+8+16+64+128",
		219 => "1+2+8+16+64+128",
		220 => "0+4+8+16+64+128",
		221 => "1+4+8+16+64+128",
		222 => "0+2+4+8+16+64+128",
		223 => "1+2+4+8+16+64+128",
		224 => "0+32+64+128",
		225 => "1+32+64+128",
		226 => "0+2+32+64+128",
		227 => "1+2+32+64+128",
		228 => "0+4+32+64+128",
		229 => "1+4+32+64+128",
		230 => "0+2+4+32+64+128",
		231 => "1+2+4+32+64+128",
		232 => "0+8+32+64+128",
		233 => "1+8+32+64+128",
		234 => "0+2+8+32+64+128",
		235 => "1+2+8+32+64+128",
		236 => "0+4+8+32+64+128",
		237 => "1+4+8+32+64+128",
		238 => "0+2+4+8+32+64+128",
		239 => "1+2+4+8+32+64+128",
		240 => "0+16+32+64+128",
		241 => "1+16+32+64+128",
		242 => "0+2+16+32+64+128",
		243 => "1+2+16+32+64+128",
		244 => "0+4+16+32+64+128",
		245 => "1+4+16+32+64+128",
		246 => "0+2+4+16+32+64+128",
		247 => "1+2+4+16+32+64+128",
		248 => "0+8+16+32+64+128",
		249 => "1+8+16+32+64+128",
		250 => "0+2+8+16+32+64+128",
		251 => "1+2+8+16+32+64+128",
		252 => "0+4+8+16+32+64+128",
		253 => "1+4+8+16+32+64+128",
		254 => "0+2+4+8+16+32+64+128",
		255 => "1+2+4+8+16+32+64+128",
		256 => "256+0",
		272 => "256+0+16",
	);
	 
	 return \%SAMFlagTableHsh;
}
########################################################################## checkFileSizeAndDefineIntervalSize
sub checkFileSizeAndDefineIntervalSize {
    
    my $fileToCheckPath = $_[0];
    my $linesToSample = $_[1];
	
	#---make sure $linesToSample is a non-zero number, if not set to 10000
	my $linesToSampleInt = int $linesToSample;
	$linesToSampleInt = 100000 if (($linesToSampleInt != $linesToSample) or ($linesToSampleInt == 0));

    #---get the filename from the path
    my @fileToCheckPathSplt = split /\//, $fileToCheckPath;
    my $fileToCheckName = $fileToCheckPathSplt[-1];
    
    print "Estimating the number of lines in $fileToCheckName.\n";
	
	#---estimate the number of lines in the file
	open (INFILE, $fileToCheckPath) || die "Can't open $fileToCheckPath.\n";
	my $tmpFilePath = $fileToCheckPath."_tmp.txt";
	system "tail -$linesToSampleInt $fileToCheckPath >$tmpFilePath";
	my $fileToCheckSize = -s "$fileToCheckPath";
   	my $tmpFileSize = -s "$tmpFilePath";
	system "rm $tmpFilePath";
   	my $fileToCheckSizeTotalLineNum = int (($fileToCheckSize/$tmpFileSize)*$linesToSampleInt);
   	print "Estimated to have ".$fileToCheckSizeTotalLineNum." lines in $fileToCheckName.\n";

	my $intervalSize = int ($fileToCheckSizeTotalLineNum/100); #---define as 
	$intervalSize = 1000000 if  ($intervalSize > 1000000);
	
	return ($fileToCheckSizeTotalLineNum, $intervalSize);

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $lineProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalLineNum = $_[3];
	my $intervalStart = $_[4];

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalLineNum - $lineProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	updateProgressBar("[$runTime]; Interval=$timeElapsed sec; EstimatedEnd=$estimatedEnd min", $lineProc, $fileTotalLineNum, 50, 5);
	#print "$lineProc lines processed. Last $intervalSize lines:".$timeElapsed." sec. Estimated end: ".$estimatedEnd." mins.\r";
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
		
}
########################################################################## updateProgressBar
sub updateProgressBar {
	
	my $strToPrint = $_[0];
	my $progressCount = $_[1];
	my $totalCount = $_[2];
	my $scaleMax = $_[3];
	my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
	my $progressPct = int (($progressCount/$totalCount)*$scaleMax);

	my $progressBar = "|";
	for my $i (1..$progressPct) {$progressBar .= ">";}
	for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
	$progressBar .= "|";

	my $extraWhiteSpaceStr = "";
	for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
	print $progressBar.$strToPrint.$extraWhiteSpaceStr."\r";

}
########################################################################## printProgressScale
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readMultiFasta
sub readMultiFasta {

	my $theRefFastaPath = $_[0];
	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	print "Reading $theRefFastaPath into a hash.\n";
	open (INFILE, $theRefFastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;
	
	return (\%fastaHsh);
}
########################################################################## outputReadLengthHitMismatch
sub outputReadLengthHitMismatch {

	my %rdLenHitMsmtchHsh = %{$_[0]};
	
	my %lengthBasedInfoHsh;
	my %allLenHsh;
	my %allNMHsh;
	my %allNHHsh;
	
	foreach my $readLength (sort {$a <=> $b} keys %rdLenHitMsmtchHsh) {
		foreach my $NM (sort {$a <=> $b} keys %{$rdLenHitMsmtchHsh{$readLength}}) {
			foreach my $NH (sort {$a <=> $b} keys %{${$rdLenHitMsmtchHsh{$readLength}}{$NM}}) {
				my $count = ${${$rdLenHitMsmtchHsh{$readLength}}{$NM}}{$NH};
				${${$lengthBasedInfoHsh{$readLength}}{"NM"}}{$NM} = 0 if (not exists ${${$lengthBasedInfoHsh{$readLength}}{"NM"}}{$NM});
				${${$lengthBasedInfoHsh{$readLength}}{"NH"}}{$NH} = 0 if (not exists ${${$lengthBasedInfoHsh{$readLength}}{"NH"}}{$NH});
				${${$lengthBasedInfoHsh{$readLength}}{"total"}}{"total"} = 0 if (not exists ${${$lengthBasedInfoHsh{$readLength}}{"total"}}{"total"});
				${${$lengthBasedInfoHsh{$readLength}}{"NH"}}{$NH} += $count;
				${${$lengthBasedInfoHsh{$readLength}}{"NM"}}{$NM} += $count;
				${${$lengthBasedInfoHsh{$readLength}}{"total"}}{"total"} += $count;
				
				$allNHHsh{$NH}++;
				$allNMHsh{$NM}++;
				$allLenHsh{$readLength}++;
			}
		}
	}
	
	open (OUTHITMSMATCH, ">$outDir/$paraTag/readLengthHitMismatch.log.txt");
	my @allLengthAry;
	foreach my $readLength (sort {$a <=> $b} keys %allLenHsh) {push @allLengthAry, $readLength;}
	my $allLengthStr = join "\t", @allLengthAry;
	print OUTHITMSMATCH "\t".$allLengthStr."\n";
	foreach my $NH (sort {$a <=> $b} keys %allNHHsh) {
		my @outLineAry;
		push @outLineAry, "NH".$NH;
		foreach my $readLength (@allLengthAry) {
			my $count = 0;
			$count = ${${$lengthBasedInfoHsh{$readLength}}{"NH"}}{$NH} if (exists ${${$lengthBasedInfoHsh{$readLength}}{"NH"}}{$NH});
			push @outLineAry, $count;
		}
		print OUTHITMSMATCH join "", ((join "\t", @outLineAry), "\n");
	}
	close OUTHITMSMATCH;
}
########################################################################## printMismatchPolyATailInfo
sub printMismatchPolyATailInfo {

	my %rdLenPolyATailCountHsh = %{$_[0]};
	my %rdLenNoPolyATailCountHsh = %{$_[1]};
	
	my %ATailLenCountHsh;
	my %nonPolyACountHsh;
	my %lengthCountHsh;
	
	my $num = keys %rdLenPolyATailCountHsh;
	if ($num > 0) {
		open (POLYATAILMSMTCH, ">$outDir/$paraTag/polyAMismatch.log.txt");
		print POLYATAILMSMTCH join "", ((join "\t", ("length", "ATailLen", "endAMismatchNum", "count")), "\n");
		foreach my $readLength (sort {$a <=> $b} keys %rdLenPolyATailCountHsh) {
			foreach my $ATailLen (sort {$a <=> $b} keys %{$rdLenPolyATailCountHsh{$readLength}}) {
				my $totalCount = 0;
				foreach my $endAMismatchNum (sort {$a <=> $b} keys %{${$rdLenPolyATailCountHsh{$readLength}}{$ATailLen}}) {
					my $count = ${${$rdLenPolyATailCountHsh{$readLength}}{$ATailLen}}{$endAMismatchNum};
					$totalCount += $count;
					print POLYATAILMSMTCH join "", ((join "\t", ($readLength, $ATailLen, $endAMismatchNum, $count)), "\n");
				}
				$ATailLenCountHsh{$ATailLen} = 0 if (not exists $ATailLenCountHsh{$ATailLen});
				$ATailLenCountHsh{$ATailLen} += $totalCount;
				my $nonPolyALen = $readLength - $ATailLen;
				$nonPolyACountHsh{$nonPolyALen} = 0 if (not exists $nonPolyACountHsh{$nonPolyALen});
				$nonPolyACountHsh{$nonPolyALen} += $totalCount;
				$lengthCountHsh{$readLength} = 0 if (not exists $lengthCountHsh{$readLength});
				$lengthCountHsh{$readLength} += $totalCount;
				
			}
		}
		
		print POLYATAILMSMTCH "\n"."ATailLen"."\t"."count"."\n";
		foreach my $ATailLen (sort {$a <=> $b} keys %ATailLenCountHsh) {
			print POLYATAILMSMTCH $ATailLen."\t".$ATailLenCountHsh{$ATailLen}."\n"; 
		}
		
		print POLYATAILMSMTCH "\n"."nonPolyALen"."\t"."count"."\n";
		foreach my $nonPolyALen (sort {$a <=> $b} keys %nonPolyACountHsh) {
			print POLYATAILMSMTCH $nonPolyALen."\t".$nonPolyACountHsh{$nonPolyALen}."\n"; 
		}

		print POLYATAILMSMTCH "\n"."polyARdFulllength"."\t"."count"."\n";
		foreach my $readLength (sort {$a <=> $b} keys %lengthCountHsh) {
			print POLYATAILMSMTCH $readLength."\t".$lengthCountHsh{$readLength}."\n"; 
		}

		print POLYATAILMSMTCH "\n"."nonPolyARdFulllength"."\t"."count"."\n";
		foreach my $readLength (sort {$a <=> $b} keys %rdLenNoPolyATailCountHsh) {
			print POLYATAILMSMTCH $readLength."\t".$rdLenNoPolyATailCountHsh{$readLength}."\n"; 
		}

		close POLYATAILMSMTCH;
	}
}
########################################################################## runSamtoolToConvertSAMToBamAndIndex
sub runSamtoolToConvertSAMToBamAndIndex {
	
	#--runSamtoolToConvertSAMToBamAndIndex($samPathHsh_ref);
	
	my %samPathHsh = %{$_[0]};
	
	#---test if fasta index if present
	my $fastaIdxPath = "$refFastaPath.fai";
	my $indexFastaCmd = "no";
	open (TEST, "$fastaIdxPath") or $indexFastaCmd = "samtools faidx $refFastaPath"; close TEST;
	
	foreach my $samPathName (keys %samPathHsh) {
		my $samPathToConvert = $samPathHsh{$samPathName};
		print "Samtools processing $samPathName\n";
		
		my $bamPath = $samPathToConvert;
		$bamPath =~ s/sam$/bam/;
		my $cmd = "samtools import $fastaIdxPath $samPathToConvert $bamPath";
		runAndCheckSerialTask("ps -ef | grep samtools | grep -v grep | grep -v perl", $cmd, $cmd);
		
		$cmd = "samtools sort $bamPath $bamPath.sorted";
		runAndCheckSerialTask("ps -ef | grep samtools | grep -v grep | grep -v perl", $cmd, $cmd);

		$cmd = "rm -f $bamPath";
		runAndCheckSerialTask("ps -ef | grep rm | grep -v grep | grep -v perl", $cmd, $cmd);

		$cmd = "samtools index $bamPath.sorted.bam";
		runAndCheckSerialTask("ps -ef | grep samtools | grep -v grep | grep -v perl", $cmd, $cmd);

		$cmd = "samtools view $bamPath.sorted.bam >$bamPath.sorted.sam";
		runAndCheckSerialTask("ps -ef | grep samtools | grep -v grep | grep -v perl", $cmd, $cmd);

		$cmd = "rm -f $samPathToConvert";
		runAndCheckSerialTask("ps -ef | grep rm | grep -v grep | grep -v perl", $cmd, $cmd);
	}
}
########################################################################## runAndCheckSerialTask
sub runAndCheckSerialTask {

	my $grepCmd = $_[0];
	my $grepStr = $_[1];
	my $cmd = $_[2];

	system (qq|$cmd &|);
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCmd`;
		sleep 1;
	}
}
########################################################################## readGff
sub readGff {

	my $gffPathToRead = $_[0];
	
	# my($nameByGeneHsh_ref, $strndByGeneHsh_ref, $cntgByGeneHsh_ref, $exonRngByGeneHsh_ref, $exonNumByCntgHsh_ref, $geneExonLenHsh_ref, $geneCtgryHsh_ref, $ctgryReadCountHsh_ref, $intronRngByGeneHsh_ref) = readGff($gffPathToRead);

	#---variables to retun
	my (%nameByGeneHsh, %strndByGeneHsh, %cntgByGeneHsh, %exonRngByGeneHsh, %exonNumByCntgHsh, %geneExonLenHsh, %geneCDSLenHsh, %geneCtgryHsh, %ctgryReadCountHsh, %CDSRngByGeneHsh);
	my (%geneByRNAHsh, %CDSCountHsh, %exonCountHsh, %geneExonLocationHsh, %intronRngByGeneHsh);

	#---read the gff
	open (INFILE, $gffPathToRead);
	print "Reading gff.\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		
		if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
			my @theLineSplt = split (/\t/, $theLine);
			my $seq = $theLineSplt[0];

			my $geneCategory = $theLineSplt[2];
			my $featureStart = $theLineSplt[3];
			my $featureEnd = $theLineSplt[4];
			my $geneStrd = $theLineSplt[6];
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my $dscrptns = $theLineSplt[8];
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				$strndByGeneHsh{$geneID} = $geneStrd;
				$cntgByGeneHsh{$geneID} = $seq;
				$nameByGeneHsh{$geneID} = $geneName;

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				# The CDS is ignored at the moment, until it reaches the point that we are looking at UTRs
				#
				#my $mRNAID = $parent;
				#my $geneID = $geneByRNAHsh{$mRNAID};
				#$CDSCountHsh{$geneID}++;
				#my $CDSCount = $CDSCountHsh{$geneID};

				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"} = $featureStart;
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} = $featureEnd;
			 	#$geneCDSLenHsh{$geneID} = 0 if $CDSCount == 1; #---define the length hashfor the 1st time
			 	#$geneCDSLenHsh{$geneID} += ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} - ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"};
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $exonID = $unqID;
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh{$RNAID};
				my $locationTag = $seq.":".$featureStart.":".$featureEnd;

				${$geneExonLocationHsh{$locationTag}}{$geneID}++;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} = $featureEnd;
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh{$RNAID} = $geneID;
				$geneCtgryHsh{$geneID} = $geneCategory;
				
				if (not(exists $ctgryReadCountHsh{$geneCategory})) {#---initialize the $geneCategory for all category 
					${$ctgryReadCountHsh{$geneCategory}}{"s"} = 0;
					${$ctgryReadCountHsh{$geneCategory}}{"a"} = 0;
				}
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;
	
	#--- check is there any ftur stored;
	my $fturNum =  keys %geneCtgryHsh;
	die "No feature stored. Likely due to conflicts in the contig names or wrong cntgFltr options\n" if ($fturNum == 0);
	
	print "Finished reading.\n";

	print "Start collapsing the alternative transcripts.\n";

	#---find the alternative transcripts that are labelled as different genes, assuming the same $seq.":".$featureStart.":".$featureEnd combination is a sign of alternative transcriptn
	my %altTrncptGeneHsh;
	foreach my $locationTag (keys %geneExonLocationHsh) {
		my $geneNum = keys %{$geneExonLocationHsh{$locationTag}};
		if ($geneNum > 1) {
			my @tmpGeneIDAry;
			foreach my $geneID (sort {$a cmp $b} keys %{$geneExonLocationHsh{$locationTag}}) {
				push @tmpGeneIDAry, $geneID;
			}
			my $geneIDString = join "\t", @tmpGeneIDAry;
			$altTrncptGeneHsh{$geneIDString}++; #---since the geneID combinations could be duplicated for many times, so have to use hash
		}
	}
	
	#---fuse the alternative transcripts as the same gene
	foreach my $geneIDString (keys %altTrncptGeneHsh) {
		my @geneIDAry = split /\t/, $geneIDString;
		my @sortedGeneIDAry = sort @geneIDAry;
		my $altTrnscptNum = @sortedGeneIDAry;
		my $newGeneID = $sortedGeneIDAry[0].".alt.".$altTrnscptNum;
		$strndByGeneHsh{$newGeneID} = $strndByGeneHsh{$sortedGeneIDAry[0]};
		$cntgByGeneHsh{$newGeneID} = $cntgByGeneHsh{$sortedGeneIDAry[0]};
		$geneCtgryHsh{$newGeneID} = $geneCtgryHsh{$sortedGeneIDAry[0]};
		$nameByGeneHsh{$newGeneID} = $nameByGeneHsh{$sortedGeneIDAry[0]};

		foreach my $geneID (@geneIDAry) {
			#----transfer the exon info to the $newGeneID and remove the old ones
			foreach my $exonID (keys %{$exonRngByGeneHsh{$geneID}}) {
				${${$exonRngByGeneHsh{$newGeneID}}{$exonID}}{"start"} = ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
				${${$exonRngByGeneHsh{$newGeneID}}{$exonID}}{"end"} = ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"};
			}
			delete $exonRngByGeneHsh{$geneID};
			delete $cntgByGeneHsh{$geneID};
			delete $strndByGeneHsh{$geneID};
		}
	}

	my $exonRngByGeneHsh_ref = \%exonRngByGeneHsh;

	my %fusedGeneHsh;
	my $fusedGeneHsh_ref = \%fusedGeneHsh;
	
	my %geneToCheckHsh = %cntgByGeneHsh;
	my $geneToCheckHsh_ref = \%geneToCheckHsh;
	
	my $overlapExon = "yes";

	while ( $overlapExon eq "yes") {
		($exonRngByGeneHsh_ref, $fusedGeneHsh_ref, $geneToCheckHsh_ref, $overlapExon) = fuseOverlappingExon($exonRngByGeneHsh_ref, $fusedGeneHsh_ref, $geneToCheckHsh_ref);
	}	
	
	my $altTrnscptGeneNum = keys %altTrncptGeneHsh;
	print "$altTrnscptGeneNum genes were found to be annotated with alternative transcripts\n";

	%exonRngByGeneHsh = %{$exonRngByGeneHsh_ref};
	my $geneNum = keys %exonRngByGeneHsh;
	print "The exon ranges of $geneNum genes have been stored.\n";

	#---store bounds, calculate the exon length
	my %boundsForIntronHsh;
	foreach my $geneID (keys %exonRngByGeneHsh) {
		my $exonCount = keys %{$exonRngByGeneHsh{$geneID}};
		my $seq = $cntgByGeneHsh{$geneID};

		${$exonNumByCntgHsh{$seq}}{$geneID} = $exonCount; #---store the ftur by contig, along side its count			

		foreach my $exonID (keys %{$exonRngByGeneHsh{$geneID}}) {
	
		 	$geneExonLenHsh{$geneID} = 0 if (not (exists $geneExonLenHsh{$geneID})); #---define the length hashfor the 1st time
			$geneExonLenHsh{$geneID} += ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} - ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
			
			if ($exonCount > 1) {
				push @{$boundsForIntronHsh{$geneID}}, ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
				push @{$boundsForIntronHsh{$geneID}}, ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"};
			}
		}
	}

	#---define the introns ranges
	my $geneIntronNum = keys %boundsForIntronHsh;

	print "$geneIntronNum gene were found to contain intron. Storing intron boundaries.\n";

	foreach my $geneID (keys %boundsForIntronHsh) {
		my @sortedBounds = sort {$a <=> $b} @{$boundsForIntronHsh{$geneID}};
		my $boundNum = @sortedBounds;
		
		my $intronNum = 0;
		for (my $i = 1; $i < ($boundNum - 1); $i += 2) {
			$intronNum++;
			my $intronID = $geneID.$intronNum;
			${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"start"} = $sortedBounds[$i] + 1;
			${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"end"} = $sortedBounds[$i+1] - 1;
		}
	}
	
	return (\%nameByGeneHsh, \%strndByGeneHsh, \%cntgByGeneHsh, \%exonRngByGeneHsh, \%exonNumByCntgHsh, \%geneExonLenHsh, \%geneCtgryHsh, \%ctgryReadCountHsh, \%intronRngByGeneHsh);
}
########################################################################## generateTwoWaysJunctionIndex
sub generateTwoWaysJunctionIndex {

	# my $twoWaysJunctIndxHsh_ref = generateTwoWaysJunctionIndex($intronRngByGeneHsh_ref, $cntgByGeneHsh_ref);
	
	my %intronRngByGeneHsh = %{$_[0]};
	my %cntgByGeneHsh = %{$_[1]};
	
	my %twoWaysJunctIndxHsh;
	
	foreach my $geneID (keys %intronRngByGeneHsh) {
		my $cntg = $cntgByGeneHsh{$geneID};
		foreach my $intronID (keys %{$intronRngByGeneHsh{$geneID}}) {
			my $intronStart = ${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"start"};
			my $intronEnd = ${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"end"};
			${${$twoWaysJunctIndxHsh{$cntg}}{"+"}}{$intronStart} = $intronEnd;
			${${$twoWaysJunctIndxHsh{$cntg}}{"-"}}{$intronEnd} = $intronStart;
		}
	}
	
	return \%twoWaysJunctIndxHsh;

}
########################################################################## getJunctionDownStreamSeq
sub checkUpAndDownStreamSeqPolyA {

	# my $polyA = checkUpAndDownStreamSeqPolyA($dynamicDownMaxAPct, $downANumRegion, $downANumMax, $downATrackRegion, $downATrackMax, $upANumRegion, $upANumMax, $cntgSeq, $readStrand, $rdName, $twoWaysJunctIndxHsh_ref, $genomicLength, $readStart, $cntg, $checkJunctionPolyA);

	my $dynamicDownMaxAPct = $_[0];
	my $downANumRegion = $_[1];
	my $downANumMax = $_[2];
	my $downATrackRegion = $_[3];
	my $downATrackMax = $_[4];
	my $upANumRegion = $_[5];
	my $upANumMax = $_[6];
	my $cntgSeq = $_[7];
	my $readStrand = $_[8];
	my $rdName = $_[9];
	my %twoWaysJunctIndxHsh = %{$_[10]};
	my $genomicLength = $_[11];
	my $readStart = $_[12];
	my $cntg = $_[13];
	my $checkJunctionPolyA = $_[14];

	my $polyA = "no";
	

	
	return $polyA;
	
	#print POLYARMREAD $theLine."\t"."UP:Z:$upRegNumCheckSeq"."\t"."DN:Z:$downRegNumCheckSeq"."\n" if ($countOnly ne "yes");#----turn off ad hoc, not print the removed read
}
########################################################################## fuseOverlappingExon
sub fuseOverlappingExon {

	my %tmpExonRngByGeneHsh = %{$_[0]};
	my %geneToCheckHsh = %{$_[1]};
	my %fusedGeneHsh = %{$_[2]};

	my (%outputExonRngByGeneHsh, $overlapExon, %outputGeneToCheckHsh);

	#----To finalize the exons, fuse the overlapping exons from alternative transcripts
	foreach my $geneID (keys %tmpExonRngByGeneHsh) {

		if (exists $geneToCheckHsh{$geneID}) {
			my %dumpExonHsh;
			foreach my $primaryExonID (keys %{$tmpExonRngByGeneHsh{$geneID}}) {

				#--- the $primaryExonID is not dumped;
				if (not (exists $dumpExonHsh{$primaryExonID})) {
					my $pStart = ${${$tmpExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"start"};
					my $pEnd = ${${$tmpExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"end"};
					my $overlapWithAllExons = "no"; 
					
					foreach my $secondaryExonID (keys %{$tmpExonRngByGeneHsh{$geneID}}) {
						
						#--- the $secondaryExonID is not dumped and not equal to primaryExonID;
						if ((not (exists $dumpExonHsh{$secondaryExonID})) and ($secondaryExonID ne $primaryExonID)){
							my $sStart = ${${$tmpExonRngByGeneHsh{$geneID}}{$secondaryExonID}}{"start"};
							my $sEnd = 	${${$tmpExonRngByGeneHsh{$geneID}}{$secondaryExonID}}{"end"};				

							#---check range overlapping
							my ($scene, $overlap) = checkRangePairOverlap($pStart, $pEnd, $sStart, $sEnd);
						
							if ($overlap eq "yes") {#--overlapping, fuse it;
							
								my @sortedRange =  sort ($pStart, $pEnd, $sStart, $sEnd);
								my $fusedStart = $sortedRange[0];
								my $fusedEnd = $sortedRange[-1];
								$dumpExonHsh{$primaryExonID}++;
								$dumpExonHsh{$secondaryExonID}++;
								$fusedGeneHsh{$geneID}++;
								my $fusedExonID = $primaryExonID."_".$secondaryExonID;
								${${$outputExonRngByGeneHsh{$geneID}}{$fusedExonID}}{"start"} = $fusedStart;
								${${$outputExonRngByGeneHsh{$geneID}}{$fusedExonID}}{"end"} = $fusedEnd;
								$overlapWithAllExons = "yes";
								print "exon $primaryExonID and $secondaryExonID were fused.\n";
							}
						}
						
						last if ($overlapWithAllExons eq "yes"); #---skip checking with the other secondary exons if the primary is fused;

					}#---end of foreach $secondaryExonID (keys %{$tmpExonRngByGeneHsh{$geneID}}) {
					
					if ($overlapWithAllExons eq "no") {#---not overlapping with other exons;
						$dumpExonHsh{$primaryExonID}++;
						${${$outputExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"start"} = ${${$tmpExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"start"};
						${${$outputExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"end"} = ${${$tmpExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"end"};
					}
					
				}#---end of if (not (exists $dumpExonHsh{$primaryExonID})) {
			} #---end of foreach $primaryExonID (keys %{$tmpExonRngByGeneHsh{$geneID}}) {

		} else {#---no need to check for overlapping, output directly;
			foreach my $exonID (keys %{$tmpExonRngByGeneHsh{$geneID}}) {
				${${$outputExonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = ${${$tmpExonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
				${${$outputExonRngByGeneHsh{$geneID}}{$exonID}}{"end"} = ${${$tmpExonRngByGeneHsh{$geneID}}{$exonID}}{"end"};
			}
		}#---end of  else {#---no need to check for overlapping, output directly;
	}#--- end of foreach my $geneID (keys %tmpExonRngByGeneHsh) {
	
	#---check if the the fused exons are still overlapping with the others;
	$overlapExon = "no";
	foreach my $geneID (keys %fusedGeneHsh) {#---check the fused genes only, for the non-fused gene, they should be possible to contain overlapping exons anymore
		foreach my $primaryExonID (keys %{$outputExonRngByGeneHsh{$geneID}}) {
			my $pStart = ${${$outputExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"start"};
			my $pEnd = ${${$outputExonRngByGeneHsh{$geneID}}{$primaryExonID}}{"end"};
			foreach my $secondaryExonID (keys %{$outputExonRngByGeneHsh{$geneID}}) {
				if ($secondaryExonID ne $primaryExonID) {
					my $sStart = ${${$outputExonRngByGeneHsh{$geneID}}{$secondaryExonID}}{"start"};
					my $sEnd = 	${${$outputExonRngByGeneHsh{$geneID}}{$secondaryExonID}}{"end"};				
					my ($scene, $overlap) = checkRangePairOverlap($pStart, $pEnd, $sStart, $sEnd);
					if ($overlap eq "yes") {$overlapExon = "yes";}
					$outputGeneToCheckHsh{$geneID}++;
				}
			}
		}
	}

	return (\%outputExonRngByGeneHsh, \%fusedGeneHsh, \%outputGeneToCheckHsh, $overlapExon);
}
########################################################################## checkRangePairOverlap
sub checkRangePairOverlap {

	my $pStart = $_[0];
	my $pEnd = $_[1];
	my $sStart = $_[2];
	my $sEnd = $_[3];

	my ($scene, $overlap);

#					The 7 scenes of overlapping and proximity 
#
#
#     case 0: complete overlapp (($sStart == $pStart) && ($sEnd == $pEnd))
#			
#     case 1: overlapHead         case 2: overlapTail	        case 3: cover		     case 4: within		case 5: prxmtyHead	     case 6: prxmtyTail
#
# primary |--------|		        |---------|	        |-------------|	              |-----|			|-----|				       	     |-------|
# secondary <=========>	         <==========>		          <=========>	            <==========>		     	    <==========>	      <==========>
#
#     ($pStart<$sStart)&&	   ($pStart>=$sStart)&&	     ($pStart<$sStart)&&      ($pStart>$sStart)&&        ($pStart<$sStart)&&		   ($pEnd>$sStart)&&
#     ($pEnd>=$sStart)&&	   ($pStart<=$sEnd)&&	     ($pEnd>$sEnd)	       ($pEnd<$sEnd)		   ($pStart<$sEnd)		   ($pStart>$sEnd)
#     ($pEnd<=$sEnd)	         ($pEnd>$sEnd)												 
#

	if  (($pStart == $sStart) && ($pEnd == $sEnd)) {$scene = 0; $overlap = "yes";}
	elsif (($pStart<=$sStart)&&($pEnd>=$sStart)&&($pEnd<=$sEnd)) {$scene = 1; $overlap = "yes";}
	elsif (($pStart>=$sStart)&&($pStart<=$sEnd)&&($pEnd>=$sEnd)) {$scene = 2; $overlap = "yes";}
	elsif (($pStart<=$sStart)&&($pEnd>=$sEnd)) {$scene = 3; $overlap = "yes";}
	elsif (($pStart>=$sStart)&&($pEnd<=$sEnd)) {$scene = 4; $overlap = "yes";}
	elsif (($pStart<$sStart)&&($pStart<$sEnd)) {$scene = 5; $overlap = "no";}
	elsif (($pEnd>$sStart)&&($pStart>$sEnd)) {$scene = 6; $overlap = "no";}
	else {#---BUG! possibly other scene?
		die "Unexpected overlapping scene. It's a bug and check your code. Program qutting.\n";
	}
	
	return ($scene, $overlap);
}
