#!/usr/bin/env perl
#2023_06_09 Purpose: reformat Bambino output to SJ format.
#Usage: perl reformat_bambino2SJ.pl -i $input -o $output -sample $sample (sample name to be used in the output file column 3, optional)
use strict;
use Getopt::Long;
use FileHandle;
use Carp;

my ($infile,$outfile,$sample);
GetOptions(
    "i=s" => \$infile,
    "o=s" => \$outfile, #optional, default $infile".sj"
    "sample=s" => \$sample, #optional 
);

#QA
croak "FATAL: you must provide infile -$infile-" unless $infile && -e $infile;
$outfile=$infile.'.sj' unless defined($outfile);

my @out;
my $fh_in=FileHandle->new($infile) or croak "cannot open $infile";
while(<$fh_in>){
	chomp;
	my $line=$_;
	next if $line=~/^NormalSample/;
	my ($NormalSample,$TumorSample,$Name,$Chr,$Pos,$Type,$Size,$Coverage,$Percent_alternative_allele,$Chr_Allele,$Alternative_Allele,$Score,$Text,$unique_alternative_ids,$reference_normal_count,$reference_tumor_count,$alternative_normal_count,$alternative_tumor_count,$count_ref_normal_fwd,$count_ref_normal_rev,$count_ref_tumor_fwd,$count_ref_tumor_rev,$count_var_normal_fwd,$count_var_normal_rev,$count_var_tumor_fwd,$count_var_tumor_rev,$alternative_fwd_count,$alternative_rev_count,$alternative_bidirectional_confirmation,$broad_coverage,$broad_reference_normal_count,$broad_reference_tumor_count,$broad_alternative_normal_count,$broad_alternative_tumor_count,$unique_alt_read_starts,$unique_alt_read_starts_fwd,$unique_alt_read_starts_rev,$avg_mapq_alternative,$somatic_or_germline,$loh_flag,$alt_to_ref_ratio_normal,$alt_to_ref_ratio_tumor,$alt_to_ref_ratio_diff_normal_tumor,$strand_skew)=split /\t/,$line;
	#QA
	for my $n ($alternative_tumor_count,$reference_tumor_count,$alternative_normal_count,$reference_normal_count){
		croak "FATAL: n must be a positive integer but not -$n-" unless ($n=~/^\d+$/ and $n >= 0);
	}
	$TumorSample=~s/\.bam$//;
	$TumorSample=~s/\_.*//;
	my $sample2=defined($sample) ? $sample : $TumorSample;
	my $Total_In_Tumor=$alternative_tumor_count+$reference_tumor_count;
	my $Total_In_Normal=$alternative_normal_count+$reference_normal_count;
	$Chr_Allele='-' unless defined($Chr_Allele);
	$Alternative_Allele='-' unless defined($Alternative_Allele);
	push(@out,(join("\t",('NA','Bambino',$sample2,$Chr,$Pos,qw/- - - - /,$alternative_tumor_count,$Total_In_Tumor,$alternative_normal_count,$Total_In_Normal,$Chr_Allele,$Alternative_Allele,$Text))));
}
$fh_in->close;

#output
my $fh_out=FileHandle->new("> $outfile") or croak "cannot write to output file $outfile";
for my $line (@out){
	print $fh_out "$line\n";
}
$fh_out->close;
