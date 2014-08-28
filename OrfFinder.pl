#! /usr/bin/perl

my $h=
"#                             *******ORF FINDER*******
# Gives info about the found ORFs according to the user definition choice.
# Including:  (only for options 1 and 2)
#				Lenght
#				Start position in the translated sequence
#				Frame
# -i|infile: Nucleotide Fasta sequence. 
# -c|choice: Custome options:  (default=1)
#		1 ORFs between Met and STOP only.
#		2 ORFs between STOP and STOP only (except for the first aac).
#		3 Translate de hole sequence to aminoacid sequence.
# -f|fasta: Create output ready for BlastP. 
#		default false, if wanted introduce the outfile desired name.
#		CAUTION: Aminoacid sequences will contain '*' in stop codon places
# -l|lenght: Stablish a treshold for lenght of ORFs, longer than 'x' (def=0)
#		For option 1 and 2, is possible to display only ORFs bigger than -l
# -b|best: Only show the longest ORF 0/1 (def=0)
#		For option 1 and 2 only.
#		Use the 6 frames. No opt -s allowed
# -s|frames: If wanted to concider both strands. 1/3/6 (def=6)
#			1 for 1rst frame
#			3 for 3 frames in the same strans
#			6 for all frames (both strands)
#
#                                             Cynthia Alexander Rascon 2010-01";

##############################ABRIR ARCHIVOS####################################
use strict;  
use warnings;
use Getopt::Long;

my %opts = ();
GetOptions (\%opts,'i|infile=s', 'c|choice=i', 'b|best=i', 'f|fasta=s',
            'l|lenght=i', 's|frames=i');


my $o=($opts{c} || 1);
my $f=($opts{f} || 0);
my $l=($opts{l} || 0);
my $b=($opts{b} || 0);
my $s=($opts{s} || 6);
open(IN,"$opts{i}") || die "Error, no input options.\n$h\n";
my @lines=<IN>;
if($f){open(OUT,">$f.fa") || die "Error\n";}

my $seq='';
my $name='';
foreach my $line (@lines){
	$line =~ s/\r|\n//g;     #para saltos de linea de windows
	chomp($line);
	if ($line=~ /^>/) { $name=$line; next;}
	if ($line=~ /^\s+/) {next;}
	$seq=$seq.$line;         #seq tendra la secuencia en una sola linea
}
$seq=uc($seq);

#######################FUNCIONES################################################

sub  Reversecomplement {         
	my $strand=$_[0];
	my $revcom= reverse ($strand);
	$revcom=~ tr/ATGC/TACG/;
	return ($revcom);                    
	
}
sub  Translation {         
	my $strand=$_[0];
	my ($i,$j,@triplets,@protein) =(0,0,(),());
	my %codons = (
	'A' => ['GCT','GCC','GCA','GCG'],
	'C' => ['TGT','TGC'],
	'D' => ['GAT','GAC'],
	'E' => ['GAA','GAG'],
	'F' => ['TTT','TTC'], 
	'G' => ['GGT','GGC','GGA','GGG'],
	'H' => ['CAT','CAC'],
	'I' => ['ATT','ATC','ATA'],
	'K' => ['AAA','AAG'],
	'L' => ['TTA','TTG','CTT','CTC','CTA','CTG'],
	'M' => ['ATG'],
	'N' => ['AAT','AAC'],
	'P' => ['CCT','CCC','CCA','CCG'],
	'Q' => ['CAA','CAG' ],
	'R' => ['CGT','CGC','CGA','CGG', 'AGA','AGG'],
	'S' => ['TCT','TCC','TCA','TCG','AGT','AGC'],
	'T' => ['ACT','ACC','ACA','ACG'],
	'V' => ['GTT','GTC','GTA','GTG'],
	'W' => ['TGG'],
	'Y' => ['TAT','TAC'],
	'*' => ['TAA','TAG','TGA']); 
	while ($strand=~ s/^(\w\w\w)//){
		$triplets[$i]=$1;
	$i++;}
	my $num= scalar (@triplets);
	for ($i=0; $i<$num;$i++){
		foreach my $aac ( keys %codons){
			foreach my $codon (@{$codons{$aac}}){
				if($codon eq $triplets[$i]){
					$protein[$j]=$aac;
					$j++;
					last;
				}}}}
	return (@protein);               
}

sub PrintProt{
	my $reftoprint=$_[0];
	my $printed='';
	my $x='';
	for $x (@{$_[0]}){        # @{...} dereferencia la referencia, refiriendose asi
                              # a el arreglo como tal
	$printed=$printed.$x}
	return ($printed);
}

####################################PROGRAMA#####################################


my $marco1=$seq;
my $marco2="NN".$seq;
my $marco3="N".$seq;
my $marco4=&Reversecomplement($seq);
my $marco5="NN".$marco4;
my $marco6="N".$marco4;
my @output1=&Translation($marco1);
my @output2=&Translation($marco2);
my @output3=&Translation($marco3);
my @output4=&Translation($marco4);
my @output5=&Translation($marco5);
my @output6=&Translation($marco6);
# @ref -- arreglo que contiene referencias a: cada uno de los arreglos que
# contienen (uno por cada marco): los aacs
my @ref=(\@output1,\@output2,\@output3,\@output4,\@output5,\@output6);
my $temp;

# c es el contador de posicion en el arreglo 'output' que tiene los tripletes
# el 1er for recorre desde la pos. 0 hasta la ultima, que es el tamaño del arreglo
# la bandera sirve para detectar la ausencia de ORF
# el segundo for parte de la posicion en que se encontro una metionina en el
#  arreglo de tripletes, hasta la ultima posicion del arreglo de tripletes (tamaño
#  - c quees la posicion actual)
my $c;
my $size;
my $flag=1;
my $pept;
my $frame=1;
my $lenght=0;
my $longest='';
my $longframe='';

if ($o==1){
	my $temp2=0;
	for my $r (@ref){   #recorre cada referencia a proteina (una por marco)
		if ($frame>$s) {print "\n\n"; last;}
		$temp2++;
		if (!$b) {print "\n\n $temp2 º reading frame:\n\t";}
        #scalar(@{$r}) da el tamaño del vector (proteina) al que apunta $r
		for ($c=0, $size=scalar(@{$r});$c<$size;$c++){
        #$r referencia a una proteina (de los posibles marcos) ->[c] es un aac
			if ($r->[$c] eq 'M'){
				for ($temp=$c, $pept='';$temp<$size;$temp++){
					$pept=$pept.$r->[$temp];
					if($r->[$temp]eq'*') {last;}
				}
				$temp=length($pept);   #para no contar 'top' en stop
				if($temp>$l){
					$flag=0;
					if ($b){ 
						if($temp>$lenght){
							$lenght=$temp;
							$longest=$pept;
							$longframe=$frame;
						}
					}
					else {print "$temp aac $pept\n\t";
                   if($f) {print OUT "$name frame $frame lenght $temp\n$pept\n";}}
				}
			}	
		}
		if ($flag){print "---No ORFs found---\n";}
		$frame++;
	}
}


if ($o==2){
	my $temp2=0;
	for my $r (@ref){   #recorre cada referencia a proteina (una por marco)
		if ($frame>$s) {print "\n\n"; last;}
		$temp2++;
		if (!$b){print "\n\n $temp2 º reading frame: \n\t";}
        #scalar(@{$r}) da el tamaño del vector (proteina) al que apunta $r
		for ($c=0, $pept='', $size=scalar(@{$r});$c<$size;$c++){   					$pept=$pept.$r->[$c];
					if($r->[$c]eq'*'||$c==($size-1)) {
						$temp=length($pept);   
						if($temp>$l){
							$flag=0;
							if ($b){ 
								if($temp>$lenght){
									$lenght=$temp;
									$longest=$pept;
									$longframe=$frame;
								}
							}
							else{print "$temp aac $pept\n\t";
							if($f) {print OUT "$name frame $frame lenght $temp\n$pept\n";}}
						}
						$pept='';
					}
		}
		if ($flag){print "---No ORFs found---\n";}
		$frame++;
	}
}

if ($b){
	print "\tLongest ORF, frame $longframe:\n\t$lenght\t$longest";
if ($f) { print OUT "$name frame $longframe lenght $lenght\n$longest\n";}}


if ($o==3){
    #es necesario hacer el paso por referencia de un arreglo a la subrutina
	$temp=PrintProt($ref[0]);
	print "\nFirst reading frame: \n\t$temp";
	if($f) {print OUT "$name _1_frame\n$temp";}
	if ($s==1) {print "\n\n"; exit;}
	$temp=PrintProt($ref[1]);
	print "\n\nSecond reading frame:\n\t$temp";
	if($f) {print OUT "\n$name _2_frame\n$temp";}
	$temp=PrintProt($ref[2]);
	print "\n\nThird reading frame:\n\t$temp";
	if($f) {print OUT "\n$name _3_frame\n$temp";}
	if ($s==3) {print "\n\n"; exit;}
	$temp=PrintProt($ref[3]);
	print "\n\nFourth reading frame:\n\t$temp";
	if($f) {print OUT "\n$name _4_frame\n$temp";}
	$temp=PrintProt($ref[4]);
	print "\n\nFifth reading frame:\n\t$temp";
	if($f) {print OUT "\n$name _5_frame\n$temp";}
	$temp=PrintProt($ref[5]);
	print "\n\nSixt reading frame:\n\t$temp";
	if($f) {print OUT "\n$name _6_frame\n$temp";}
}
	
print "\n\n";
close (IN);
if ($f) {close (OUT);}