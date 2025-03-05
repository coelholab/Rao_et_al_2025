#! /usr/bin/perl


use strict;
use warnings;
use Readonly;
use Getopt::Std;
use File::Basename;

Readonly::Scalar my $B_TRUE  => (1 == 1);
Readonly::Scalar my $B_FALSE => (1 == 0);
Readonly::Scalar my $S_EMPTY => '';

sub trim
{
  my @localcpy = @_;
  my @trimcpy = ();
  foreach (@localcpy)
  {
    next if( ! defined );
    s/^\s+//; s/\s+$//;
    push @trimcpy , $_ if( $_ ne $S_EMPTY );
  }
  return (wantarray)?@trimcpy:$trimcpy[0];
}

sub usage
{
  print STDERR <<_EOM_
 ---------------------------------------------------------------------
 Usage: $^X $0 > -f <fasta file>
 ---------------------------------------------------------------------
_EOM_
;
  exit 1;
}

my %opts;
getopts( 'f:e:a:h', \%opts );

my $clones_file;


&usage if(defined $opts{h});

$clones_file = $opts{f} if(defined $opts{f} && -f $opts{f});

&usage if(! defined $clones_file);

my @fasta;


my $header;
my %sequences;
open FH, "$clones_file" or die "Error: reading FASTA file [$!]";
while(my $line = <FH>)
{
	chomp $line;
	
	if($line =~ /^[>]/){             #split the fasta file based on the '>' identifier at the beginning of each FASTA annotation
		$line = &trim($line);
		$header = $line;
		$header =~ s/^[>]//g;
		#printf "%s\n",$header;
		#printf "%s\n",$header;
	}
	
	if(!($line =~ /^[>]/)){             #split the fasta file based on the '>' identifier at the beginning of each FASTA annotation
		my @residues = split(//, $line);
		$sequences{$header} = [] if(! defined $sequences{$header});
		push @{$sequences{$header}}, @residues;
	}
	
	
}
close FH;

foreach my $h (keys %sequences){
	
	#printf "%s\n",$h;
	
}


my %mutations;

if(defined $sequences{'Germline'}){
	
	#printf "Yes\n";
	
}



for(my $a = 0; $a < scalar(@{$sequences{'Germline'}}); $a++){
	
	my $pos = $a + 1;
	
	foreach my $seq (keys %sequences){
	
		if($seq ne '>Germline'){
			
			if($sequences{$seq}[$a] ne $sequences{'Germline'}[$a]){
				
				my $mut = $sequences{'Germline'}[$a].$pos.$sequences{$seq}[$a];
				
				
									
				$mutations{$mut} = [] if(! defined $mutations{$mut});
				push @{$mutations{$mut}}, $seq;
					
					
				
				
				

				
			}
			
		}
	
	
	}
	
	
}


foreach my $m (keys %mutations){
	
	printf "%s\t",$m;
	
	foreach my $s (@{$mutations{$m}}){
		
		printf "%s,",$s;
		
		
	}
	
	printf "\n";
	
	
}



