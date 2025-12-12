#!/usr/bin/perl -w

#  COPYRIGHT NOTICE
#  
#  Silico - A Perl Molecular Toolkit
#  Copyright (C) 2008-25 David K. Chalmers and Benjamin P. Roberts,
#  Department of Medicinal Chemistry, Monash University
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#<
#! silico_definitions.pm
#? Define Silico global constants.
#. This module contains no subroutines but defines severable global arrays.
#  These are:
#.
#. @Atomic_Elements - an array containing the periodic table.
#. @Atomic_Masses - array containing atomic masses.
#. @Amino_Acids - the standard amino acids.
#. @Amino_Acids - the standard amino acids.
#. @Amino_Acids1 - the standard amino acids
#. @Nucleic_Acids - the standard nucleic acids.
#. @Nucleic_Acids1 - the standard nucleic acids - lower case one letter code.
#. %Amino_Nucleic_Acids - hash with AA and NA three-letter code as key with one-letter codes
#  as value.
#. %Amino_Nucleic_Acids1 - hash with AA and NA one-letter code as key with three-letter codes
#  as value.
#. @Vdw_Radii - VdW radii for determinining if two atoms are bonded
#. $Revision: 1.32.2.3.2.1 $
#>

use strict;
package Silico;

fundamental_constants();
atomic_elements();
silico_vdw_radii();
amino_acid_and_nucleotide_base_codes();

sub fundamental_constants {
	
	#<
	#? Set various mathematical and physical fundamental constants
	#  as read-only values (so they can't be reassigned by a wayward
	#  use of "=" instead of, for example, "==").
	#; Requires: nothing
	#; Returns: nothing (sets package global arrays and hashes)
	#>
	
	use vars qw($Avogadro $Elementary_Charge $Pi $R $Tet);
	
	# Avogadro's Number
	*Silico::Avogadro = \6.02214199e23;
	
	# Elementary Charge
	*Silico::Elementary_Charge = \1.60217646e-19;
	
	# Pi
	*Silico::Pi = \3.14159265358979;
	
	# Gas Constant
	*Silico::R = \8.31447215;
	
	# Tetrahedral angle
	*Silico::Tet = \109.471220634;
	
	return undef;
}


sub atomic_elements {
	
	#<
	#? Prepare an array atomic element symbols, with the offset corresponding to
	#  the atomic number (dummy atom is 0).  Also prepare a hash with the element
	#  name as the key and the element number as the value
	#; Requires: nothing
	#; Returns: nothing (sets package global arrays and hashes)
	#>
	
	return if $Silico::Atomic_Elements->[0];

	#	1   2   3   f1  f2  f3  f4  f5  f6  f7  f8  f9  f10 f11 f12 f13 f14 4   5   6   7   8   9   10  11  12  13  14  15  16  17  18
	@Silico::Atomic_Elements = qw(
	
		Du
	
		H                                                                                                                           He
		Li  Be                                                                                                  B   C   N   O   F   Ne
		Na  Mg                                                                                                  Al  Si  P   S   Cl  Ar
		K   Ca  Sc                                                          Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr
		Rb  Sr  Y                                                           Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe
		Cs  Ba  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn
		Fr  Ra  Ac  Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr	
	);
	
	my $i = 0;
	foreach (@Silico::Atomic_Elements) {
		$Silico::Atomic_Elements{$_} = $i;
		++$i;
	}
	
}

sub atomic_masses {
	
	#<
	#? Return a list of atomic masses, with the offset corresponding to
	#  the atomic number (dummy atom is 0).
	#. Atomic masses are taken from the CRC handbook, Standard Atomic
	#  Weights, 1981.
	#; Requires: nothing 
	#; Returns: nothing (sets package global arrays and hashes)

	use vars qw(@Atomic_Masses);
	
	# Group:   1          2          3         f1         f2         f3         f4         f5         f6         f7         f8         f9         f10        f11        f12        f13        f14         4          5          6          7          8          9          10         11         12         13         14         15         16         17         18
	
	@Silico::Atomic_Masses = qw(
				  
		  0
		
		  1.00794                                                                                                                                                                                                                                                                                                                                              4.00260
		  6.941      9.01218                                                                                                                                                                                                                                                                           10.81      12.011     14.0067    15.9994    18.998403  20.179
		 22.98977   24.305                                                                                                                                                                                                                                                                             26.98154   28.0855    30.97376   32.06      35.453     39.948
		 39.0983    40.08      44.9559                                                                                                                                                              47.88      50.9415    51.996     54.938     55.847     58.9332    58.69      63.456     65.38      69.72      72.59      74.9216    78.96      79.904     83.80
		 85.4678    87.62      88.9059                                                                                                                                                              91.22      92.9064    95.94      98        101.07     102.9055   106.42     107.8682   112.41     114.82     118.69     121.75     127.60     126.9045   131.29
		132.9054   137.33     138.9055   140.12     140.9077   144.24     145        150.36     151.96     157.25     158.9254   162.50     164.9304   167.26     168.9342   173.04     174.967    178.49     180.9479   183.85     186.207    190.2      192.22     195.08     196.9665   200.59     204.383    207.2      208.9804   209        210        222
		223        226.0254   227.0278   232.0381   321.0359   238.0289   237.0482   244        243        247        247        251        252        257        258        259        260
	

	);
	
	return undef;
}

sub amino_acid_and_nucleotide_base_codes {
	
	#<
	#? Prepare lists of three-letter and one-letter codes for
	#  amino acids and nucleotide bases. Also prepare hashes
	#  mapping the three-letter and one-letter codes to each
	#  other.
	#; Requires: nothing
	#; Returns: nothing  (sets package global arrays and hashes)
	#>
	
	return if $Silico::Amino_Acids->[0];
	
	
	# Note: Where there are multiple defnitions, the primary definition must come last in this list (eg HIS, CYS, XXX)
	@Silico::Amino_Acids  = qw(--- ACE ALA ARG ASP ASN CYX CYM CYS GLU GLN GLY HID HIE HIP HSD HSE HSP HIS ILE LEU LYS MET NME NMA SER THR PHE PRO TRP TYR VAL XXX);
	@Silico::Amino_Acids1 = qw(-   X   A   R   D   N   C   C   C   E   Q   G   H   H   H   H   H   H   H   I   L   K   M   X   X   S   T   F   P   W   Y   V   X  );
	
	# Note: PDB residue names 'C G T A U' refer to ribosenucleic acids (RNA) 'DC DG DT DA DU' refer to 
	# 2'-deoxyribosenucleic acis (DNA)

	@Silico::Nucleic_Acids  = qw(--- DC DG DT DA DU C G T A U XXX);
	@Silico::Nucleic_Acids1 = qw(-   c  g  t  a  u  c g t a u X);
	
	# %Amino_Nucleic_Acids (three letter keys and 1-letter values) and %Amino_Nucleic_Acids1
	# (one letter keys and three letter values) are hashes of all usual
	# Nucleic and Amino Acids that can be used to identify usual residues and
	# to translate one letter code to three letter code and vice-versa
	
	my $i = 0;
	foreach (@Silico::Amino_Acids)  {
		
		$Silico::Amino_Nucleic_Acids{$_} = $Silico::Amino_Acids1[$i];
		$Silico::Amino_Nucleic_Acids1{$Silico::Amino_Acids1[$i]} = $_;
		
		++$i;
	}

	$i = 0;
	foreach (@Silico::Nucleic_Acids) {
		
		$Silico::Amino_Nucleic_Acids{$_} = $Silico::Nucleic_Acids1[$i];
		$Silico::Amino_Nucleic_Acids1{$Silico::Nucleic_Acids1[$i]} = $_;
	
		++$i;
	}
	
	return undef;
}

sub silico_vdw_radii {
	

	#<
	#? Create a list of Silico van der Waals radii, one for each atom.
	#. This data should be used only to determine whether or not to
	#  create a bond when generating connectivity information.
	#. Some elements have their radii artifically set to zero, so as to
	#  prevent bonds from being created to these elements.
	#. The hydrogen atom's radius is reduced to 0.2 Angstroms.
	#. Bonds to (between?) C, H, N, O, P and S are determined elsewhere
	#  with reference to the expected lengths for bonds to these atoms.
	#; Requires: nothing
	#; Returns: nothing
	#>
	
	use vars qw(@Vdw_Radii);

	return if $Silico::Vdw_Radii[0];
	
	#	1    2    3    f1   f2   f3   f4   f5   f6   f7   f8   f9   f10  f11  f12  f13  f14  4    5    6    7    8    9    10   11   12   13   14   15   16   17   18
	@Silico::Vdw_Radii = qw(
		
		0.00
		
		0.20                                                                                                                                                       0.00
		0.00 0.00                                                                                                                         2.08 1.65 1.54 1.45 1.35 0.00
		0.00 0.00                                                                                                                         2.05 2.00 1.90 2.00 1.81 0.00
		0.00 0.00 1.70                                                                       2.55 2.55 2.55 2.55 2.55 2.55 2.55 2.55 2.55 1.70 1.70 2.00 2.00 2.10 0.00
        	0.00 0.00 1.70                                                                       1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 2.20 2.20 2.50 0.00
		0.00 0.00 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 1.70 0.00
		0.00 0.00 1.70 1.70 1.70 1.70 1.70 0.90 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 1.20
	
	);
	
	return undef;
}

sub standard_bond_lengths {
	
	#<
	#? Create a hash containing standard chemical bond lengths.
	#  The hash key contains the bond: two atoms, between which
	#  is a character denoting the bond order.
	#  Dash ('-'): single bond
	#  Equals sign ('='): double bond
	#  Percent sign ("%"): triple bond
	#. At the moment, bond lengths are given twice: one for A-B, the
	#  other for B-A. This is because we can't be sure in advance
	#  which order the atoms will be given in.
	#. Bond lengths were sourced from
	#  http://chemviz.ncsa.uiuc.edu/content/doc-resources-bond.pdf,
	#  January 2007.
	#. Note: The bond length for Sb-Cl is for antimony trichloride,
	#  not antimony pentachloride.
	#; Requires: nothing
	#; Returns: nothing
	#>
	
	use vars qw(%Bond_Lengths);

	return if $Silico::Bond_Lengths->[0];
	
	# Group:			1		13		14		15		16		17
	%Silico::Bond_Lengths = qw(	H-H	0.740	B-H	1.190	C-H	1.090	N-H	1.010	O-H	0.960	F-H	0.920
					H-B	1.190
					H-C	1.090			C-C	1.540	N-C	1.470	O-C	1.430	F-C	1.350
									C=C	1.340	N=C	1.290	O=C	1.200
									C%C	1.200	N%C	1.160	O%C	1.130
					H-N	1.010			C-N	1.470	N-N	1.450	O-N	1.400	F-N	1.360
									C=N	1.290	N=N	1.250	O=N	1.210
									C%N	1.160	N%N	1.100
					H-O	0.960			C-O	1.430	N-O	1.400	O-O	1.480	F-O	1.420
									C=O	1.200	N=O	1.210	O=O	1.210
									C%O	1.130
					H-F	0.920			C-F	1.350	N-F	1.360	O-F	1.420	F-F	1.420
					H-Si	1.480			C-Si	1.850			O-Si	1.630	F-Si	1.600
					H-P	1.440			C-P	1.840			O-P	1.630	F-P	1.540
													O=P	1.500
					H-S	1.340			C-S	1.820					F-S	1.560
									C=S	1.600			O=S	1.430
					H-Cl	1.270	B-Cl	1.750	C-Cl	1.770	N-Cl	1.750
					H-Ge	1.530			C-Ge	1.950					F-Ge	1.680
					H-As	1.520							O-As	1.780	F-As	1.710
					H-Se	1.460
					H-Br	1.410			C-Br	1.940
					H-Sn	1.700			C-Sn	2.160
					H-Te	1.700
					H-I	1.610			C-I	2.140					F-I	1.910
									C-Pb	2.300
									Si-H	1.480	P-H	1.440	S-H	1.340	Cl-H	1.270
															Cl-B	1.750
									Si-C	1.850	P-C	1.840	S-C	1.820	Cl-C	1.770
													S=C	1.600
															Cl-N	1.750
									Si-O	1.630	P-O	1.630
											P=O	1.500	S=O	1.430
									Si-F	1.600	P-F	1.540	S-F	1.560
									Si-Si	2.330			S-Si	2.000	Cl-Si	2.020
											P-P	2.210			Cl-P	2.030
													S=P	1.860
									Si-S	2.000			S-S	2.050	Cl-S	1.820
											P=S	1.860	S=S	1.490
									Si-Cl	2.020	P-Cl	2.030	S-Cl	2.070	Cl-Cl	1.990
															Cl-Ge	2.100
															Cl-As	2.160
									Si-Br	2.150
															Cl-Sn	2.330
															Cl-Sb	2.320
									Si-I	2.430	P-I	1.840			Cl-I	2.320
															Cl-Pb	2.420
									Ge-H	1.530	As-H	1.520	Se-H	1.460	Br-H	1.410
									Ge-C	1.950					Br-C	1.940
											As-O	1.780
									Ge-F	1.680	As-F	1.710
															Br-Si	2.150
									Ge-Cl	2.100	As-Cl	2.160
									Ge-Ge	2.410					Br-!Ge	2.300
											As-As	2.430			Br-As	2.330
													Se=Se	2.150
									Ge-Br	2.300	As-Br	2.330			Br-Br	2.280
															Br-Sn	2.500
											As-I	2.540
									Sn-H	1.700			Te-H	1.700	I-H	1.610
									Sn-C	2.160					I-C	2.140
															I-F	1.910
									Sn-Cl	2.330	Sb-Cl	2.320			I-Cl	2.320
															I-As	2.540
									Sn-Br	2.500
															I-Sn	2.700
									Sn-I	2.700					I-I	2.670
															I-Pb	2.790
									Pb-C	2.300
									Pb-Cl	2.420
									Pb-I	2.790
		);
		
	return undef;
}

sub maximum_valences {
	
	#<
	#? Returns a hash containing, for each element, the maximum usual neutral
	#  valence for that element.
	#; Requires: nothing
	#; Returns: hash
	#>
	
	my $maxvalences;
	
	%$maxvalences = qw(	Du 0
				H  1 He 0 Li 1 Be 2 B  4 C  4 N  4 O  2 F  1 Ne 0
				Na 1 Mg 2 Al 4 Si 4 P  5 S  6 Cl 1 Ar 0	K  1 Ca 2
				Sc 3 Ti 4 V  5 Cr 6 Mn 7 Fe 3 Co 3 Ni 3 Cu 2 Zn 2
				Ga 3 Ge 4 As 5 Se 6 Br 1 Kr 0 Rb 1 Sr 2 Y  3 Zr 4
				Nb 5 Mo 6 Tc 7 Ru 8 Rh 4 Pd 4 Ag 1 Cd 2 In 3 Sn 4
				Sb 5 Te 6 I  1 Xe 0 Cs 1 Ba 2 La 3 Ce 4 Pr 3 Nd 3
				Pm 3 Sm 3 Eu 3 Gd 3 Tb 4 Dy 3 Ho 3 Er 3 Tm 3 Yb 3
				Lu 3 Hf 4 Ta 5 W  6 Re 7 Os 8 Ir 6 Pt 6 Au 7 Hg 2
				Tl 3 Pb 4 Bi 5 Po 4 At 1 Rn 0 Fr 1 Ra 2 Ac 3 Th 4
				Pa 5 U  6);
	
	return $maxvalences;
}

sub acceptable_neutral_valences {
	
	#<
	#? Determines, for each element, a list of acceptable neutral
	#  valences. The ANV is like an oxidation state, but it
	#  doesn't take into account atom electronegativities.
	#; Requires: nothing (creation of package global variable)
	#; Returns: nothing (creation of package global variable)
	#>
	
	use vars qw($Valences);

	# We don't need to do this again if it has been done
	return if ($Silico::Valences->{H});
	
	@{$Silico::Valences->{H}}  = (  1);
	@{$Silico::Valences->{He}} = (0);
	@{$Silico::Valences->{Li}} = (  1);
	@{$Silico::Valences->{Be}} = (    2);
	@{$Silico::Valences->{B}}  = (      3);
	@{$Silico::Valences->{C}}  = (        4);
	@{$Silico::Valences->{N}}  = (      3,  5);
	@{$Silico::Valences->{O}}  = (    2);
	@{$Silico::Valences->{F}}  = (  1);
	@{$Silico::Valences->{Ne}} = (0);
	@{$Silico::Valences->{Na}} = (  1);
	@{$Silico::Valences->{Mg}} = (    2);
	@{$Silico::Valences->{Al}} = (      3);
	@{$Silico::Valences->{Si}} = (        4);
	@{$Silico::Valences->{P}}  = (      3,  5);
	@{$Silico::Valences->{S}}  = (    2,  4,  6);
	@{$Silico::Valences->{Cl}} = (  1,  3,  5,  7);
	@{$Silico::Valences->{Ar}} = (0);
	@{$Silico::Valences->{K}}  = (  1);
	@{$Silico::Valences->{Ca}} = (    2);
	@{$Silico::Valences->{Sc}} = (      3);
	@{$Silico::Valences->{Ti}} = (    2,  4);
	@{$Silico::Valences->{V}}  = (    2,3,4,5);
	@{$Silico::Valences->{Cr}} = (      3,    6);
	@{$Silico::Valences->{Mn}} = (    2,3,4,  6,7);
	@{$Silico::Valences->{Fe}} = (    2,3);
	@{$Silico::Valences->{Co}} = (    2,3);
	@{$Silico::Valences->{Ni}} = (    2,3);
	@{$Silico::Valences->{Cu}} = (  1,2);
	@{$Silico::Valences->{Zn}} = (    2);
	@{$Silico::Valences->{Ga}} = (      3);
	@{$Silico::Valences->{Ge}} = (        4);
	@{$Silico::Valences->{As}} = (      3,  5);
	@{$Silico::Valences->{Se}} = (    2,  4,  6);
	@{$Silico::Valences->{Br}} = (  1,  3,  5,  7);
	@{$Silico::Valences->{Kr}} = (0);
	@{$Silico::Valences->{Rb}} = (  1);
	@{$Silico::Valences->{Sr}} = (    2);
	@{$Silico::Valences->{Y}}  = (      3);
	@{$Silico::Valences->{Zr}} = (        4);
	@{$Silico::Valences->{Nb}} = (      3,  5);
	@{$Silico::Valences->{Mb}} = (    2,3,4,5,6);
	@{$Silico::Valences->{Tc}} = (  	    7);
	@{$Silico::Valences->{Ru}} = (    2,3,4,  6,  8);
	@{$Silico::Valences->{Rh}} = (    2,3,4);
	@{$Silico::Valences->{Pd}} = (    2,  4);
	@{$Silico::Valences->{Ag}} = (  1);
	@{$Silico::Valences->{Cd}} = (    2);
	@{$Silico::Valences->{In}} = (      3);
	@{$Silico::Valences->{Sn}} = (    2,  4);
	@{$Silico::Valences->{Sb}} = (      3,  5);
	@{$Silico::Valences->{Te}} = (    2,  4,  6);
	@{$Silico::Valences->{I}}  = (  1,  3,  5,  7);
	@{$Silico::Valences->{Xe}} = (0);
	@{$Silico::Valences->{Cs}} = (  1);
	@{$Silico::Valences->{Ba}} = (    2);
	@{$Silico::Valences->{La}} = (      3);
	@{$Silico::Valences->{Ce}} = (      3,4);
	@{$Silico::Valences->{Pr}} = (      3);
	@{$Silico::Valences->{Nd}} = (      3);
	@{$Silico::Valences->{Pm}} = (      3);
	@{$Silico::Valences->{Sm}} = (      3);
	@{$Silico::Valences->{Eu}} = (    2,3);
	@{$Silico::Valences->{Gd}} = (      3);
	@{$Silico::Valences->{Tb}} = (      3,4);
	@{$Silico::Valences->{Dy}} = (      3);
	@{$Silico::Valences->{Ho}} = (      3);
	@{$Silico::Valences->{Er}} = (      3);
	@{$Silico::Valences->{Tm}} = (      3);
	@{$Silico::Valences->{Yb}} = (    2,3);
	@{$Silico::Valences->{Lu}} = (      3);
	@{$Silico::Valences->{Hf}} = (        4);
	@{$Silico::Valences->{Ta}} = (      3,4,5);
	@{$Silico::Valences->{W}}  = (  1,2,3,4,5,6);
	@{$Silico::Valences->{Re}} = (  1,2,3,4,5,6,7);
	@{$Silico::Valences->{Os}} = (  1,2,3,4,5,6,7,8);
	@{$Silico::Valences->{Ir}} = (    2,3,4,  6);
	@{$Silico::Valences->{Pt}} = (  1,2,3,4,5,6);
	@{$Silico::Valences->{Au}} = (  1,2,3,4,5,6,7);
	@{$Silico::Valences->{Hg}} = (  1,2);
	@{$Silico::Valences->{Tl}} = (  1,  3);
	@{$Silico::Valences->{Pb}} = (    2,  4);
	@{$Silico::Valences->{Bi}} = (      3,  5);
	@{$Silico::Valences->{Po}} = (    2,   4);
	@{$Silico::Valences->{At}} = (  1,  3,  5,  7);
	@{$Silico::Valences->{Rn}} = (0);
	@{$Silico::Valences->{Fr}} = (  1);
	@{$Silico::Valences->{Ra}} = (    2);
	@{$Silico::Valences->{Ac}} = (      3);
	@{$Silico::Valences->{Th}} = (        4);
	@{$Silico::Valences->{Pa}} = (    2,3,4,5);
	@{$Silico::Valences->{U}}  = (      3,4,5,6);
	
	return undef;
}

sub electronegativities {
	
	#<
	#? Create an array of Pauling electronegativities. Atoms are given
	#  in order of atomic number (dummy atom is 0).
	#; Requires: nothing 
	#; Returns: nothing (creation of package global array)
	#>
	
	use vars qw(@Electronegativities);

	# We don't need to do this again if it has been done
	return if ($Silico::Electronegativities->[0]);
	
	#	 1    2    3   f1   f2   f3   f4   f5   f6   f7   f8   f9   f10  f11  f12  f13  f14   4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
	
	@Silico::Electronegativities = qw(
	
		0.00
		
		2.20                                                                                                                                                       0.00
		0.98 1.57                                                                                                                         2.04 2.55 3.04 3.44 3.98 0.00
		0.93 1.31                                                                                                                         1.61 1.90 2.19 2.58 3.16 0.00
		0.82 1.00 1.36                                                                       1.54 1.63 1.66 1.55 1.83 1.88 1.91 1.90 1.65 1.81 2.01 2.18 2.55 2.96 3.00
        	0.82 0.95 1.22                                                                       1.33 1.60 2.16 1.90 2.20 2.28 2.20 1.93 1.69 1.78 1.96 2.05 2.10 2.66 2.60
		0.79 0.89 1.10 1.12 1.13 1.14 1.13 1.17 1.20 1.20 1.20 1.22 1.23 1.24 1.25 1.10 1.27 1.30 1.50 2.36 1.90 2.20 2.20 2.28 2.54 2.00 1.62 2.33 2.02 2.00 2.20 0.00
		0.70 0.90 1.10 1.30 1.50 1.38 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 1.20
	
	);
	
	return undef;
}

sub space_groups {
	
	#<
	#? Generate a list of space groups to provide crystal type information.
	#  This is broadly based on the information provided in the Tripos Mol2
	#  CRYSIN record.
	#. Space groups are given in the format AAAB.
	#  AAA is a three-digit code corresponding to the crystal type number
	#  in the International Tables for X-Ray Crystallography. For example,
	#  type 1 is given as "001", while type 230 is given as "230".
	#  B is a one-digit code corresponding to the particular unit cell
	#  setting (orientation of the axes, etc.), which can be used to
	#  determine the full space group.
	#  Conventions, based on Sybyl (which is believed to be based on the
	#  aforementioned International Tables):
	#. Groups 001 and 002: only one possible setting (thus, 0011 and 0021).
	#. Groups 003 to 015: Up to four possible settings:
	#  - Unique b-axis: settings a'  b'  c' (1) and  c'  b'  a' (2)
	#  - Unique c-axis: settings a'  c' -b' (3) and -b'  c'  a' (4)
	#. Groups 016 to 074: Up to six possible settings:
	#  a  b  c (1)
	#  c  a  b (2)
	#  b  c  a (3)
	#  a -c  b (4)
	#  b  a -c (5)
	# -c  b  a (6)
	#
	#. Groups 075 to 142: Two possible settings:
	#     a     b  c (1)
	#  (a+b) (b-a) c (2)
	#  (a-b) (b+a) c (2)
	#
	#. Groups 143 to 230: Only one possible setting
	#
	#. Further conventions are used for text identification (Silico only):
	# a dash, "-", means that the preceding character has a macron over it
	# the text "sub" indicates that the following character is subscripted
	#; Requires: nothing 
	#; Returns: nothing (creation of package global array)
	
	use vars qw(%Space_Groups);

	# We don't need to do this again if it has been done
	return if ($Silico::SpaceGroups->[0]);
	
	%Silico::Space_Groups = qw(	0011	P1
					0021	P1-
					0031	P2					0033	P2
					0041	P2sub1					0043	P2sub1
					0051	C2		0052	A2		0053	B2		0054	A2
					0061	Pm					0063	Pm
					0071	Pc		0072	Pa		0073	Pb		0074	Pa
					0081	Cm		0082	Am		0083	Bm		0084	Am
					0091	Cc		0092	Aa		0093	Bb		0094	Aa
					0101	P2/m					0103	P2/m
					0111	P2sub1/m				0113	P2sub1/m
					0121	C2/m		0122	A2/m		0123	B2/m		0124	A2/m
					0131	P2/c		0132	P2/a		0133	P2/b		0134	P2/a
					0141	P2sub1/c	0142	P2sub1/a	0143	P2sub1/b	0144	P2sub1/a
					0151	C2/c		0152	A2/a		0153	B2/b		0154	A2/a
					0161	P222
					0171	P222sub1	0172	P2sub122	0173	P22sub12
					0181	P2sub12sub12	0182	P22sub12sub1	0183	P2sub122sub1
					0191	P2sub12sub12sub1
					0201	C222sub1	0202	A2sub122	0203	B22sub12
					0211	C222		0212	A222		0213	B222
					0221	F222
					0231	I222
					0241	I2sub12sub12sub1
					0251	Pmm2		0252	P2mm		0253	Pm2m
					0261	Pmc2sub1	0262	P2sub1ma	0263	Pb2sub1m	0264	Pm2sub1b	0265	Pcm2sub1	0266	P2sub1am
					0271	Pcc2		0272	P2aa		0273	Pb2b
					0281	Pma2		0282	P2mb		0283	Pc2m		0284	Pm2a		0285	Pbm2		0286	P2cm
					0291	Pca2sub1	0292	P2sub1ab	0293	Pc2sub1b	0294	Pb2sub1a	0295	Pbc2sub1	0296	P2sub1ca
					0301	Pnc2		0302	P2na		0303	Pb2n		0304	Pn2b		0305	Pcn2		0306	P2an
					0311	Pmn2sub1	0312	P2sub1mn	0313	Pn2sub1m	0314	Pm2sub1n	0315	Pnm2sub1	0316	P2sub1nm
					0321	Pba2		0322	P2cb		0323	Pc2a
					0331	Pna2sub1	0332	P2sub1nb	0333	Pc2sub1n	0334	Pn2sub1a	0335	Pbn2sub1	0336	P2sub1cn
					0341	Pnn2		0342	P2nn		0343	Pn2n
					0351	Cmm2		0352	A2mm		0353	Bm2m
					0361	Cmc2sub1	0362	A2sub1ma	0363	Bb2sub1m	0364	Bm2sub1b	0365	Ccm2sub1	0366	A2sub1am
					0371	Ccc2		0372	A2aa		0373	Bb2b
					0381	Amm2		0382	B2mm		0383	Cm2m		0384	Am2m		0385	Bmm2		0386	C2mm
					0391	Abm2		0392	B2cm		0393	Cm2a		0394	Ac2m		0395	Bma2		0396	C2mb
					0401	Ama2		0402	B2mb		0403	Cc2m		0404	Am2a		0405	Bbm2		0406	C2cm
					0411	Aba2		0412	B2cb		0413	Cc2a		0414	Ac2a		0415	Bba2		0416	C2cb
					0421	Fmm2		0422	F2mm		0423	Fm2m
					0431	Fdd2		0432	F2dd		0433	Fd2d
					0441	Imm2		0442	I2mm		0443	Im2m
					0451	Iba2		0452	I2cb		0453	Ic2a
					0461	Ima2		0462	I2mb		0463	Ic2m		0464	Im2a		0465	Ibm2		0466	I2cm
					0471	Pmmm
					0481	Pnnn
					0491	Pccm		0492	Pmaa		0493	Pbmb
					0501	Pban		0502	Pncb		0503	Pcna
					0511	Pmma		0512	Pbmm		0513	Pmcm		0514	Pmam		0515	Pmmb		0516	Pcmm
					0521	Pnna		0522	Pbnn		0523	Pncn		0524	Pnan		0525	Pnnb		0526	Pcnn
					0531	Pmna		0532	Pbmn		0533	Pncm		0534	Pman		0535	Pnmb		0536	Pcnm
					0541	Pcca		0542	Pbaa		0543	Pbcb		0544	Pbab		0545	Pccb		0546	Pcaa
					0551	Pbam		0552	Pmcb		0553	Pcma
					0561	Pccn		0562	Pnaa		0563	Pbnb
					0571	Pbcm		0572	Pmca		0573	Pbma		0574	Pcmb		0575	Pcam		0576	Pmab
					0581	Pnnm		0582	Pmnn		0583	Pnmn
					0591	Pmmn		0592	Pnmm		0593	Pmnm
					0601	Pbcn		0602	Pnca		0603	Pbna		0604	Pcnb		0605	Pcan		0606	Pnab
					0611	Pbca					0613	Pcab
					0621	Pnma		0622	Pbnm		0623	Pmcn		0624	Pnam		0625	Pmnb		0626	Pcmn
					0631	Cmcm		0632	Amma		0633	Bbmm		0634	Bmmb		0635	Ccmm		0636	Amam
					0641	Cmca		0642	Abma		0643	Bbcm		0644	Bmab		0645	Ccmb		0646	Acam
					0651	Cmmm		0652	Ammm		0653	Bmmm
					0661	Cccm		0662	Amaa		0663	Bbmb
					0671	Cmma		0672	Abmm		0673	Bmcm		0674	Bmam		0675	Cmmb		0676	Acmm
					0681	Ccca		0682	Abaa		0683	Bbcb		0684	Bbab		0685	Cccb		0686	Acaa
					0691	Fmmm
					0701	Fddd
					0711	Immm
					0721	Ibam		0722	Imcb		0723	Icma
					0731	Ibca					0733	Icab
					0741	Imma		0742	Ibmm		0743	Imcm		0744	Imam		0745	Immb		0746	Icmm
					0751	P4		0752	C4
					0761	P4sub1		0762	C4sub1
					0771	P4sub2		0772	C4sub2
					0781	P4sub3		0782	C4sub3
					0791	I4		0792	F4
					0801	I4sub1		0802	F4sub1
					0811	P4-		0812	C4-
					0821	I4-		0822	F4-
					0831	P4/m		0832	C4/m
					0841	P4sub2/m	0842	C4sub2/m
					0851	P4/n		0852	C4/a
					0861	P4sub2/n	0862	C4sub2/a
					0871	I4/m		0872	F4/m
					0881	I4sub1/a	0882	F4sub1/d
					0891	P422		0892	C422
					0901	P42sub12	0902	C422sub1
					0911	P4sub122	0912	C4sub122
					0921	P4sub12sub12	0922	C4sub122sub1
					0931	P4sub222	0932	C4sub222
					0941	P4sub22sub12	0942	C4sub222sub1
					0951	P4sub322	0952	C4sub322
					0961	P4sub32sub12	0962	C4sub322sub1
					0971	I422		0972	F422
					0981	I4sub122	0982	F4sub122
					0991	P4mm		0992	C4mm
					1001	P4bm		1002	C4mb
					1011	P4sub2cm	1012	C4sub2mc
					1021	P4sub2nm	1022	C4sub2mn
					1031	P4cc		1032	C4cc
					1041	P4nc		1042	C4cn
					1051	P4sub2mc	1052	C4sub2cm
					1061	P4sub2bc	1062	C4sub2cb
					1071	I4mm		1072	F4mm
					1081	I4cm		1082	F4mc
					1091	I4sub1md	1092	F4sub1dm
					1101	I4sub1cd	1102	F4sub1dc
					1111	P4-2m		1112	C4-m2
					1121	P4-2c		1122	C4-c2
					1131	P4-2sub1m	1132	C4-m2sub1
					1141	P4-2sub1c	1142	C4-c2sub1
					1151	P4-m2		1152	C4-2m
					1161	P4-c2		1162	C4-2c
					1171	P4-b2		1172	C4-2b
					1181	P4-n2		1182	C4-2b
					1191	I4-m2		1192	F4-2m
					1201	I4-c2		1202	F4-2c
					1211	I4-2m		1212	F4-m2
					1221	I4-2d		1222	F4-d2
					1231	P4/mmm		1232	C4/mmm
					1241	P4/mcc		1242	C4/mcc
					1251	P4/nbm		1252	C4/amb
					1261	P4/nnc		1262	C4/acn
					1271	P4/mbm		1272	C4/mmb
					1281	P4/mnc		1282	C4/mcn
					1291	P4/nmm		1292	C4/amm
					1301	P4/ncc		1302	C4/acc
					1311	P4sub2/mmc	1312	C4sub2/mcm
					1321	P4sub2/mcm	1322	C4sub2/mmc
					1331	P4sub2/nbc	1332	C4sub2/acb
					1341	P4sub2/nnm	1342	C4sub2/amn
					1351	P4sub2/mbc	1352	C4sub2/mbc
					1361	P4sub2/mnm	1362	C4sub2/mmn
					1371	P4sub2/nmc	1372	C4sub2/acm
					1381	P4sub2/ncm	1382	C4sub2/amc
					1391	I4/mmm		1392	F4/mmm
					1401	I4/mcm		1402	F4/mmc
					1411	I4sub1/amd	1412	F4sub1/ddm
					1421	I4sub1/acd	1422	F4sub1/ddc
					1431	P3
					1441	P3sub1
					1451	P3sub2
					1461	R3
					1471	P3-
					1481	R3-
					1491	P3sub12
					1501	P3sub21
					1511	P3sub112
					1521	P3sub121
					1531	P3sub212
					1541	P3sub221
					1551	R32
					1561	P3m1
					1571	P31m
					1581	P3c1
					1591	P31c
					1601	R3m
					1611	R3c
					1621	P3-1m
					1631	P3-1c
					1641	P3-m1
					1651	P3-c1
					1661	R3-m
					1671	R3-c
					1681	P6
					1691	P6sub1
					1701	P6sub5
					1711	P6sub2
					1721	P6sub4
					1731	P6sub3
					1741	P6-
					1751	P6/m
					1761	P6sub3/m
					1771	P622
					1781	P6sub122
					1791	P6sub522
					1801	P6sub222
					1811	P6sub422
					1821	P6sub322
					1831	P6mm
					1841	P6cc
					1851	P6sub3cm
					1861	P6sub3mc
					1871	P6-m2
					1881	P6-c2
					1891	P6-2m
					1901	P6-2c
					1911	P6/mmm
					1921	P6/mcc
					1931	P6sub3/mcm
					1941	P6sub3/mmc
					1951	P23
					1961	F23
					1971	I23
					1981	P2sub13
					1991	I2sub13
					2001	Pm3
					2011	Pn3
					2021	Fm3
					2031	Fd3
					2041	Im3
					2051	Pa3
					2061	Ia3
					2071	P432
					2081	P4sub232
					2091	F432
					2101	F4sub132
					2111	I432
					2121	P4sub332
					2131	P4sub132
					2141	I4sub132
					2151	P4-3m
					2161	F4-3m
					2171	I4-3m
					2181	P4-3n
					2191	F4-3c
					2201	I4-3d
					2211	Pm3m
					2221	Pn3n
					2231	Pm3n
					2241	Pn3m
					2251	Fm3m
					2261	Fm3c
					2271	Fd3m
					2281	Fd3c
					2291	Im3m
					2301	Ia3-d
					2311	P2sub1/n
	);
	
	return undef;
}

sub primes {

	#<
	#? Create an array of the first thousand prime numbers: @$Silico::primes
	#; Requires: nothing 
	#; Returns: nothing (creation of package global array)
	#>

	@$Silico::primes = qw(
              2      3      5      7     11     13     17     19     23     29 
	     31     37     41     43     47     53     59     61     67     71 
	     73     79     83     89     97    101    103    107    109    113 
	    127    131    137    139    149    151    157    163    167    173 
	    179    181    191    193    197    199    211    223    227    229 
	    233    239    241    251    257    263    269    271    277    281 
	    283    293    307    311    313    317    331    337    347    349 
	    353    359    367    373    379    383    389    397    401    409 
	    419    421    431    433    439    443    449    457    461    463 
	    467    479    487    491    499    503    509    521    523    541 
	    547    557    563    569    571    577    587    593    599    601 
	    607    613    617    619    631    641    643    647    653    659 
	    661    673    677    683    691    701    709    719    727    733 
	    739    743    751    757    761    769    773    787    797    809 
	    811    821    823    827    829    839    853    857    859    863 
	    877    881    883    887    907    911    919    929    937    941 
	    947    953    967    971    977    983    991    997   1009   1013 
	   1019   1021   1031   1033   1039   1049   1051   1061   1063   1069 
	   1087   1091   1093   1097   1103   1109   1117   1123   1129   1151 
	   1153   1163   1171   1181   1187   1193   1201   1213   1217   1223 
	   1229   1231   1237   1249   1259   1277   1279   1283   1289   1291 
	   1297   1301   1303   1307   1319   1321   1327   1361   1367   1373 
	   1381   1399   1409   1423   1427   1429   1433   1439   1447   1451 
	   1453   1459   1471   1481   1483   1487   1489   1493   1499   1511 
	   1523   1531   1543   1549   1553   1559   1567   1571   1579   1583 
	   1597   1601   1607   1609   1613   1619   1621   1627   1637   1657 
	   1663   1667   1669   1693   1697   1699   1709   1721   1723   1733 
	   1741   1747   1753   1759   1777   1783   1787   1789   1801   1811 
	   1823   1831   1847   1861   1867   1871   1873   1877   1879   1889 
	   1901   1907   1913   1931   1933   1949   1951   1973   1979   1987 
	   1993   1997   1999   2003   2011   2017   2027   2029   2039   2053 
	   2063   2069   2081   2083   2087   2089   2099   2111   2113   2129 
	   2131   2137   2141   2143   2153   2161   2179   2203   2207   2213 
	   2221   2237   2239   2243   2251   2267   2269   2273   2281   2287 
	   2293   2297   2309   2311   2333   2339   2341   2347   2351   2357 
	   2371   2377   2381   2383   2389   2393   2399   2411   2417   2423 
	   2437   2441   2447   2459   2467   2473   2477   2503   2521   2531 
	   2539   2543   2549   2551   2557   2579   2591   2593   2609   2617 
	   2621   2633   2647   2657   2659   2663   2671   2677   2683   2687 
	   2689   2693   2699   2707   2711   2713   2719   2729   2731   2741 
	   2749   2753   2767   2777   2789   2791   2797   2801   2803   2819 
	   2833   2837   2843   2851   2857   2861   2879   2887   2897   2903 
	   2909   2917   2927   2939   2953   2957   2963   2969   2971   2999 
	   3001   3011   3019   3023   3037   3041   3049   3061   3067   3079 
	   3083   3089   3109   3119   3121   3137   3163   3167   3169   3181 
	   3187   3191   3203   3209   3217   3221   3229   3251   3253   3257 
	   3259   3271   3299   3301   3307   3313   3319   3323   3329   3331 
	   3343   3347   3359   3361   3371   3373   3389   3391   3407   3413 
	   3433   3449   3457   3461   3463   3467   3469   3491   3499   3511 
	   3517   3527   3529   3533   3539   3541   3547   3557   3559   3571 
	   3581   3583   3593   3607   3613   3617   3623   3631   3637   3643 
	   3659   3671   3673   3677   3691   3697   3701   3709   3719   3727 
	   3733   3739   3761   3767   3769   3779   3793   3797   3803   3821 
	   3823   3833   3847   3851   3853   3863   3877   3881   3889   3907 
	   3911   3917   3919   3923   3929   3931   3943   3947   3967   3989 
	   4001   4003   4007   4013   4019   4021   4027   4049   4051   4057 
	   4073   4079   4091   4093   4099   4111   4127   4129   4133   4139 
	   4153   4157   4159   4177   4201   4211   4217   4219   4229   4231 
	   4241   4243   4253   4259   4261   4271   4273   4283   4289   4297 
	   4327   4337   4339   4349   4357   4363   4373   4391   4397   4409 
	   4421   4423   4441   4447   4451   4457   4463   4481   4483   4493 
	   4507   4513   4517   4519   4523   4547   4549   4561   4567   4583 
	   4591   4597   4603   4621   4637   4639   4643   4649   4651   4657 
	   4663   4673   4679   4691   4703   4721   4723   4729   4733   4751 
	   4759   4783   4787   4789   4793   4799   4801   4813   4817   4831 
	   4861   4871   4877   4889   4903   4909   4919   4931   4933   4937 
	   4943   4951   4957   4967   4969   4973   4987   4993   4999   5003 
	   5009   5011   5021   5023   5039   5051   5059   5077   5081   5087 
	   5099   5101   5107   5113   5119   5147   5153   5167   5171   5179 
	   5189   5197   5209   5227   5231   5233   5237   5261   5273   5279 
	   5281   5297   5303   5309   5323   5333   5347   5351   5381   5387 
	   5393   5399   5407   5413   5417   5419   5431   5437   5441   5443 
	   5449   5471   5477   5479   5483   5501   5503   5507   5519   5521 
	   5527   5531   5557   5563   5569   5573   5581   5591   5623   5639 
	   5641   5647   5651   5653   5657   5659   5669   5683   5689   5693 
	   5701   5711   5717   5737   5741   5743   5749   5779   5783   5791 
	   5801   5807   5813   5821   5827   5839   5843   5849   5851   5857 
	   5861   5867   5869   5879   5881   5897   5903   5923   5927   5939 
	   5953   5981   5987   6007   6011   6029   6037   6043   6047   6053 
	   6067   6073   6079   6089   6091   6101   6113   6121   6131   6133 
	   6143   6151   6163   6173   6197   6199   6203   6211   6217   6221 
	   6229   6247   6257   6263   6269   6271   6277   6287   6299   6301 
	   6311   6317   6323   6329   6337   6343   6353   6359   6361   6367 
	   6373   6379   6389   6397   6421   6427   6449   6451   6469   6473 
	   6481   6491   6521   6529   6547   6551   6553   6563   6569   6571 
	   6577   6581   6599   6607   6619   6637   6653   6659   6661   6673 
	   6679   6689   6691   6701   6703   6709   6719   6733   6737   6761 
	   6763   6779   6781   6791   6793   6803   6823   6827   6829   6833 
	   6841   6857   6863   6869   6871   6883   6899   6907   6911   6917 
	   6947   6949   6959   6961   6967   6971   6977   6983   6991   6997 
	   7001   7013   7019   7027   7039   7043   7057   7069   7079   7103 
	   7109   7121   7127   7129   7151   7159   7177   7187   7193   7207
	   7211   7213   7219   7229   7237   7243   7247   7253   7283   7297 
	   7307   7309   7321   7331   7333   7349   7351   7369   7393   7411 
	   7417   7433   7451   7457   7459   7477   7481   7487   7489   7499 
	   7507   7517   7523   7529   7537   7541   7547   7549   7559   7561 
	   7573   7577   7583   7589   7591   7603   7607   7621   7639   7643 
	   7649   7669   7673   7681   7687   7691   7699   7703   7717   7723 
	   7727   7741   7753   7757   7759   7789   7793   7817   7823   7829 
	   7841   7853   7867   7873   7877   7879   7883   7901   7907   7919 
	);
	
	return undef;
	
}
			
return 1;
