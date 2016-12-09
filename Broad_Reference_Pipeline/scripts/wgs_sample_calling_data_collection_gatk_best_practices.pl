#!/usr/bin/perl 
###########################################################
# Data Collection Script for GATK Best Practices Pipeline #
###########################################################

if (scalar(@ARGV) < 5) {
    die("Usage: SampleName NumThreads InputDataDirectory TempOutputDirectory profiling \n[if profiling is enabled, then the following is required]: collectstatspath interval stats \n");
}
my $sample =  $ARGV[0];
my $numThreads = $ARGV[1];
my $inDataDir = $ARGV[2];
my $tmpDir = "".$ARGV[3];

my $profiling = $ARGV[4]; #by default profiling is turned ON if invoked from the workflow profiler

# arguments for collect_stats
my $collectstatspath = $ARGV[5];
my $interval = $ARGV[6]; # by default sampling interval is 30s from the workflow profiler.
my $stats = $ARGV[7];

#my $numLanes =$ARGV[1];
my $called = "$0 @ARGV";

my $numLanes = 0;
my $sampleprefix = $sample.'_'.$numThreads.'T';

# INPUT FASTQ FILES
#OTHER FORMATS FOR FQ: "_1.fastq.gz; #"_1.fastq";
my $fqFile1, $fqFile2;
$fqFile1 = $inDataDir.$sample."_1.fq";
$fqFile2 = $inDataDir.$sample."_2.fq";

# Pipeline executables and its directories
### ENTER THE CORRECT PATH TO THE FOLLOWING 3 VARIABLES ###
my $broadDir = '/mnt/app_hdd/scratch/aprabh2/Broad_Ref_Pipeline/broad';
my $homosapiensrefgenomeDir = '/mnt/app_hdd/scratch/aprabh2/Broad_Ref_Pipeline/homosapiensrefgenome';
# TOOLS
my $bwaDir = "$broadDir/bwa";
my $bwa = "$bwaDir/bwa";
my $gatkDir = "$broadDir/gatk-protected/target";
my $gatk = "$gatkDir/GenomeAnalysisTK.jar";
my $picardDir ="$broadDir/picard/dist";
my $picard = "$picardDir/picard.jar";
# HOMOSAPIENSREFGENOME
my $refgenomeFastaFile = "$homosapiensrefgenomeDir/Homo_sapiens_assembly19.fasta";
my $refgenomeBwtFile = "$homosapiensrefgenomeDir/Homo_sapiens_assembly19.fasta.bwt";
my $dbSNPvcf = "$homosapiensrefgenomeDir/Homo_sapiens_assembly19.dbsnp.vcf";
my $dbSNPindel = "$homosapiensrefgenomeDir/Homo_sapiens_assembly19.known_indels.vcf";
# EXOME TARGET INTERVALS
my $exome_targets_intervals = "$homosapiensrefgenomeDir/nexterarapidcapture_exome_uniqueintervals.bed";
unless(-d $inDataDir) {
    die("Error: The InputDataDirectory $inDataDir doesn't exist\n");
}

unless(-d $tmpDir) {
    die("Error: The TempOutputDirectory $tmpDir doesn't exist\n");
}

# Output file names for each stage of the pipeline
my $baseName = $sample;
my $baseNameLane = $baseName.'_'.$numLanes.'L_'.$numThreads.'T';
my $bwamem_samFile = $tmpDir.$baseNameLane.".sam";
my $sort_bamFile = $tmpDir.$baseName."_sorted.bam";
my $duplicateMetricsFile = $tmpDir.$baseName."_dup.metrics";
my $bamDupRemFile = $tmpDir.$baseName."_dupRem.bam";
my $bamRealignFile = $tmpDir.$baseName."_realign.bam";
my $realnInterval = $tmpDir.$baseName."_realn.intervals";
my $finalBam = $tmpDir.$baseName."_final.bam";
my $HCvcf = $tmpDir.$baseName."_HaplotypeCaller.vcf";
my $genomeImportFile = $refgenomeFastaFile.".fai";
my $readGroupHeader = "\@RG\\tID:$baseNameLane\\tLB:$baseName\\tSM:$baseName\\tPL:ILLUMINA";
my $recalOut = $tmpDir .$baseName."_recal.grp";

my $dryRun = 0;

my $pwd = `pwd`;
chomp $pwd;
my $host = `hostname`;
chomp $host;
my $uname = `whoami`;
chomp $uname;
my $runningTime = time;
my $commandsfile = $tmpDir.$uname."_".$sampleprefix."_processing.log";
my $outputfile = $tmpDir.$uname."_".$sampleprefix."_output.log";
open(LOG,">$commandsfile");
print LOG "#$called (version $version) in $pwd on $host.\n";
print LOG "#Started at ".`date +"%F %T"`."\n";
print LOG "#temporary files created in $tmpDir\n";
my $procTime = time;
my $procFlag = 0;
sub run_and_log {
    my $command = $_[0];
    my $execute = !$dryRun;
    my $exitValue = 0; #a command we don't run is considered successful
    my $redirect;

    if (@_ >1){
	$execute=!$_[1];
    }

#several of the programs like to output to STDERR, so we link that to log file
    $redirect = "1>>$outputfile 2>&1";
#that is because we redirect STDOUT in many cases, so let's not mess with it.
    $redirect = "2>>$outputfile" if $command =~ m/>/;
#    $redirect = "" if $command =~ m/>/;

    $command = $command." ".$redirect;

    if ($procFlag == 0){
	$procFlag++;
    } else {
	$procTime = time - $procTime;
	printf LOG "#Processing Time %02d:%02d:%02d\n",int($procTime /3600),int(($procTime % 3600) /60),int($procTime %60);
	$procTime = time;
    }

    print LOG "#not run\n" if !$execute;
    print LOG "#".`date +"%F %T"`;
    print LOG $command."\n";

    $exitValue = system($command) if $execute;

#necessary if we use `` instead of system()
#$exitValue = $? >>8;

##If the command failed, we want to stop it here.
    if ($exitValue != 0){
	my $error = "Command failed with return value $exitValue : $command \n";
	print LOG $error;
	close LOG;
	die $error;
    }
}

#############
# Profiling #
#############
sub Start_profiling {
my ($tag) = @_;
if ($profiling) { 
system("$collectstatspath $stats -d $interval -td $tmpDir -n $sampleprefix -tag $tag -l 5 -u 1 -s 600 &");
}
}

sub Stop_Profiling {
if ($profiling) { 
system("$collectstatspath --kill-all");
}
}

##########################################
# STAGES OF GATK BEST PRACTICES PIPELINE # 
##########################################
my $stage_tag=BwaMem;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "$bwa mem -t $numThreads -Ma -R \'$readGroupHeader\' $refgenomeFastaFile $fqFile1 $fqFile2 > $bwamem_samFile";
Stop_Profiling();
sleep(60);

my $stage_tag=SortSam;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -jar $picard SortSam I=$bwamem_samFile O=$sort_bamFile SO=coordinate CREATE_INDEX=true";
Stop_Profiling();
sleep(60);

my $stage_tag=MarkDuplicates;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -jar $picard MarkDuplicates I=$sort_bamFile O=$bamDupRemFile M=$duplicateMetricsFile CREATE_INDEX=true TMP_DIR=$tmpDir";
Stop_Profiling();
sleep(60);

my $stage_tag=RealignerTargetCreator;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -jar $gatk -T RealignerTargetCreator -nt $numThreads -R $refgenomeFastaFile -o $realnInterval -known:indels,vcf $dbSNPindel -I $bamDupRemFile";
Stop_Profiling();
sleep(60);

my $stage_tag=IndelRealigner;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $gatk -T IndelRealigner -R $refgenomeFastaFile -targetIntervals $realnInterval -known:indels,vcf $dbSNPindel -I $bamDupRemFile -o $bamRealignFile --filter_bases_not_stored";
Stop_Profiling();
sleep(60);

my $stage_tag=BaseRecalibrator;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx4g -jar $gatk -T BaseRecalibrator -I $bamRealignFile -R $refgenomeFastaFile -knownSites:mask,vcf $dbSNPvcf -o $recalOut";
Stop_Profiling();
sleep(60);

my $stage_tag=PrintReads;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -jar $gatk -T PrintReads -R $refgenomeFastaFile -I $bamRealignFile -BQSR $recalOut -o $finalBam";  
Stop_Profiling();
sleep(60);

my $stage_tag=HaplotypeCaller;
Start_profiling($stage_tag);
print "$stage_tag\n";
run_and_log "java -Xmx8g -jar $gatk -T HaplotypeCaller -R $refgenomeFastaFile -I $finalBam -o $HCvcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000";
Stop_Profiling();
sleep(60);

$runningTime = time - $runningTime;

printf LOG "#done in %02d:%02d:%02d\n",int($runningTime /3600),int(($runningTime % 3600) /60),int($runningTime %60);

exit 0;
