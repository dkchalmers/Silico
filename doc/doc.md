## Directory: /vcp1/people/david/bin/Silico1.15/bin 


## bilayer_builder 


Build an approximate lipid bilayer using periodic boundary conditions. 
Molecules can be emedded in the bilayer. 


For MD packages that require anisotropic simulations to have the bilayer in 
the XY plane, use the -xy flag. 


By default, the bilayer is built in the XZ plane, using instances of a 
supplied lipid file. 


Lipid molecules are reoriented along their polar axis, unless the -no flag is 
specified 


Bilayers can be built with multiple constituent species. To do this, use 
multiple passes. That is, instances of one species are built into a bilayer, 
and this bilayer is then "embedded" into a bilayer of the second species, and 
so on. 


The distance between the two halves of the bilayer is controlled by the -s 
flag. This should be set to approximately the length of the lipid in the Y 
dimension. 


Usage: bilayer_builder <lipid_file> -x <cellx> -y <celly> -z <cellz> -n1 
<number of lipids layer1> [ -n2 <number of lipids layer2> -e <molecule to 
embed> ] 


Flags: 


-e &lt;filename&gt; Embed molecule filename <br>
-n1 &lt;number&gt; Number of instances of lipid to use to create monolayer 1 (top) <br>
-n2 &lt;number&gt; Number of instances of lipid to use to create monolayer 2 (bottom) <br>
-x &lt;number&gt; | <br>
-y &lt;number&gt; | Dimensions of the cell in each direction <br>
-z &lt;number&gt; | <br>
-xy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Build bilayer in the xy plane <br>
-no&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not reorient the polar axis of embedded lipid moleculesq <br>
-trans&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate the embedded molecule to the cell centre <br>
-rand &lt;number&gt; Number of lipid bonds to randomise <br>
-es &lt;number&gt; Additional distance addedd between the mid-points of the bilayers <br>
-si &lt;factor&gt; Initial scale factor for atom-atom packing distance (default 0.8) <br>
-sm &lt;factor&gt; Minimum scale factor for atom-atom packing distance (default 0.4) <br>
-incr &lt;number&gt; Number of unsuccesful moleucle packing trials before reducing the scale factor. Default 100 <br>
-ignh &lt;0 or 1&gt; Ignore hydrogen atoms (default 1) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 and later 


Revision: 1.27.2.1.2.19 


## checkbox 


Test to see if molecule(s) is (are) inside the box defined by box.pdb 


Use: checkbox <ligandfile> -b<boxfile> 


Used for programs like DOCK which do not always dock a molecule INSIDE the 
specified box. The default (Dock format) box is 'box.pdb'. 'Good' structures 
are written to <filebase>_ibox>.fmt and structures with atoms outside the box 
are written to <filebase>_obox.fmt 


Flags: 


-ih&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore hydrogen atoms <br>
-b&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; box location (default box.pdb) <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force overwriting of output files <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.17.2.1.2.2 


## chop_cell 


Chop down simulation cell by deleting atoms outside a newly defined cell 


Cell centred on protein atoms and then any molecule with atoms projecting 
beyond the new cell dimensions is deleted 


Flags: 


-x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; X dimension of cell used for chopping <br>
-y&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Y dimension of cell used for chopping <br>
-z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Z dimension of cell used for chopping <br>
-xr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Reduction of x dimension for cell output <br>
-yr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Reduction of x dimension for cell output <br>
-zr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Reduction of x dimension for cell output <br>
-hex&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create hexagonal cell for gromacs .gro output (not for pdb format) <br>
-ip&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore protein molecules (those with an identifiable amino acid CA atom) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.1 


Created: DKC November 1999 and later. 


## Cluster 


Cluster molecules using the very simple 'cluster' subroutine. 


Molecules are assumed to be in order of increasing energy. All molecules 
within the RMS threshold of molecule 1 are placed in the first cluster. The 
first molecule to fall outside that set is used to make a second cluster and 
so on. By default only the first member of each cluster is retained. 


Files are clustered individually unless the -a flag is set. 


Molecules are renamed to include cluster number 


Use: cluster <file> [<flags>] 


Writes out a representative structure from each cluster or all molecules 
sorted into clusters if the -r flag is specified 


Flags: 


-t &lt;number&gt; Clustering threshold (default: 2) <br>
-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Treat all input files as one set <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Retain all molecules <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 


Revision: 1.22.2.1.2.4 


## dist_filter 


Filter docking output based on protein-ligand contact with a specified proten 
residue, (residue atom, ligand atom element) 


Use options -i and -l <ligand_name> for output containing both ligand and 
protein in the same structure (e.g. Glide induced fit output) 


Flags: 


-asp &lt;specifier&gt; Atom specifier for protein atom(s)s <br>
-asl &lt;specifier&gt; Atom specifire for ligand atom(s) <br>
-p &lt;protfile&gt; Protein file name <br>
-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Input structure contains both protein and ligand (e.g. From Glide IFD docking) <br>
-max &lt;dist&gt; Filter out ligands with dist greater than this value. <br>
-min &lt;dist&gt; Filter out ligands with dist less than this value. <br>
-s &lt;string&gt; Skip ligands with names starting with 'string'  <br>

#### Atom specifier examples: 


ANAME:CA&nbsp;&nbsp;&nbsp;&nbsp; All atoms called 'CA'. <br>
ANAME:CA,CB,CG Returns all atoms called 'CA' or 'CB' or 'CD'. <br>
ANAME:CA,RESNAME:TRP All atoms called 'CA' in all residues called TRP <br>
ANAME:CA,SUBID:4 All atoms called 'CA' in residue number 4 <br>
ELEMENT:!H All nonhydrogen atoms <br>
SEGID:PROT,SEGID:LIG All atoms with the SEGID set to PROT or LIG <br>

Successive atom specifications can be made. Each is separated by a '|' 


ANAME:CB|ANAME:CA|ANAME:CD Returns all atoms called 'CB', 'CA', 'CD'. <br>
ANAME:CA,RESNUM:1|ANAME:CA,RESNUM:4 Returns 'CA' atoms from residue 1 and 4. <br>

Atom specifiers are case sensitive. 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2 


Created: DKC 2022 


## dock_sort 


Cluster and sort Glide, ICM, Vina or Gold docking results according to their 
fitness value (Glide gscoree, ICM Score or Gold fitness or liagnd efficiency) 


Processing steps: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Read docking output  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Discard molecules that have low tautomer probablility as determined by ligprep. The threshold is determined by the -tp flag  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Simplify the molecule name (for Gold output)  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Calculate the SMILES string and uses this to sort atoms into a canonical order  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Identify multiple poses of a the same molecule using a canonical SMILES string (unless the -name flag is used, then the molecule name is used)  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Cluster molecule using an RMS threshold (default 1 Angstrom)  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sort molecule atoms using a SMILES string to improve the superimposition of molecules that come from different sources and have different atom ordering  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Keep only the best pose of each ligand (if -best option is used)  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Calculate ligand efficiency  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sort molecules by docking score  <br>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Calculate enrichment curve and, if decoys have been used, a ROC (receiver operating characteristic) curve. The -decoy option can be used to specify a string used in the compound name to identify decoy compounds  <br>



Use: 


dock_sort glide_pv.maegz -o sdf -best 


Sort docking output files by energy. Options to retain only the single best 
pose of each compound (identified by SMILES string). 


Processes output from Glide, Gold, ICM and Autodock Vina. Autodock Vina output 
should be converted to .sdf format first. Using pdbqt files directly results 
in loss of atomuc element information. 


dock_sort glide_pv.maegz -o sdf -best -name 


As above but instead compounds are identified by name. 


Flags: 


### Scoring 


-f &lt;string&gt; String containing fitness value (eg glide_score or Gold.Goldscore.Fitness). Standard Glide and Gold scoring is handled automatically <br>
-ligeff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort by ligand efficiency <br>
-effm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Efficiency calculation method. le = score/natoms, r2 = score/sqrt(natoms)RAng <br>

### Retaining molecules 


-best&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Retain only single best scoring pose of each compound <br>
-t &lt;value&gt; Rms cluster threshold (default 1 Ang) <br>
-smin &lt;value&gt; Minimum value for fitness score <br>
-smax &lt;value&gt; Maximum value for fitness score <br>
-lemin &lt;value&gt; Minimum value for score/HA <br>
-lemax &lt;value&gt; Maximum value for score/HA <br>
-mwmin &lt;value&gt; Minimum MW <br>
-mwmin &lt;value&gt; Maximum MW <br>
-pct &lt;value&gt; Percentage of molecules to retain after ranking <br>
-n &lt;value&gt; Number of molecules to retain after ranking <br>
-nd &lt;value&gt; Number of molecules to retain after ranking in 'divide' directory <br>
-tp &lt;value&gt; Minimum acceptable probability for tautomers from Ligprep <br>

### Compound identification 


-name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Identify molecules by molecule name (otherwise SMILES is used) <br>
-nosmiles&nbsp;&nbsp; Speed things up by skipping SMILES generation (uses identify-by-name) <br>
-cmax &lt;value&gt; Maximum cutoff value (useful for Gold Fitness) <br>
-cmin &lt;value&gt; Minimum cutoff value (useful for Glide Score) <br>
-tp &lt;value&gt; Minimum acceptable probability for tautomers from ligprep <br>

### Contacts and H-bonds 


-c &lt;string1_string2_...&gt; Only keep ligands that contact these atoms. String format chain:resnum:aname[:functional-group_atom-element] Use quotes around multiple strings. Functional groups are defined by subroutine mol_label_functional_group(). E.g., AROMATIC_C, CARBONYL_O, S_AMIDE_N. Use the program 'mol_label_fg' to find functional groups present in a molecules. <br>
-cd &lt;value&gt; Contact atom distance cutoff (Angstroms). Default 4 Ang <br>
-ci&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Invert contact so that ligands contacting atoms are discarded <br>
-hb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Identify protein-ligand hydrogen bonds <br>
-hbmin &lt;num&gt; Filter output based on minimum &lt;num&gt; of protein-ligand hydrogen bonds <br>
-prot &lt;filename&gt; Protein file (if not available from docking output) <br>

### Compound identification 


-name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Identify molecules by molecule name. Otherwise, compounds are identified by SMILES. Note that SMILES identification currently does not encode stereochemistry <br>

### Ignoring compounds 


-ignore &lt;string&gt; Ignore compounds that start with this name. They will not be included as decoys or actives. Multiple compound names can be included by separating with spaces and surrounding with quotes. Eg -ignore "alpha beta gamma" will ignore all names starting with alpha, beta or gamma <br>
-ignores &lt;"smiles_smiles_..."&gt; Ignore compounds with specified SMILES. They will not be included as decoys or actives. Multiple compounds can be included as for -ignore <br>

### Decoys and Actives 




** Not tested lately. Probably does not work ** 




-decoy &lt;string&gt; String to identify decoy compounds. Should be present in the names of decoy compounds <br>
-afile &lt;name&gt; File containing list of active compounds, one per line <br>
-nactive &lt;int&gt; Total number of actives in the docked library <br>
-ndecoy &lt;int&gt; Total number of decoys in the docked library <br>
-act &lt;field&gt; Calculate correlation coefficents for this activity field in SDF data <br>

### Other analysis 


-hist&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make histogram files for score and ligeff values <br>

### Output 


-obase &lt;base&gt; Base name for output files and directories <br>
-divide&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write individual files to a new output directory <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Overwrite existing output directory for -divide option <br>
-term &lt;format&gt; Gnuplot output format [postscript, png] (default postscript) <br>
-rename &lt;field&gt; Rename output molecules using SDF field <br>

Revision: 1.26.2.5.2.28 


## ens_randomise_order 


Randomise the order of molecules in an ensemble 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.2.2.1.2.1 


Created: DKC November 1999 and later. 


## ens_sort 


Sort a set of molecules by some property 


-mw&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort molecules by molecular weight <br>
-ha&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort molecules by heavy atoms <br>
-s &lt;field&gt; Sort molecules using SDF data field <br>

Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.3.2.1.2.2 


Created: DKC November 1999 and later. 


## extract_lig_prot 


Split pdb file into protein and non-covalently bound ligands 


Structures are split up using connectivity 


Molecules with fewer than 'maxatom' and more 'minatom' atoms are assumed to be 
ligands (This approach is a little simplistic. It assumes that there are no 
breaks in the protein chain. However it IS able to extract peptide ligands 
from proteins) 


Flags: 


-maxatoms&nbsp;&nbsp; (Maximum ligand size) <br>
-minatoms&nbsp;&nbsp; (Minimum ligand size) <br>

Created: DKC 2000 


Revision: 1.15-80-g67aad2a 


## file_rename 


Rename files using a perl regular expression substitution 


For example file_rename -a _new <filename> will replace the string '_new' with 
'' and the file clozapine_new.mol2 becomes clozapine.mol2 


Use: 


file_rename <inputfiles> -a <exp1> -b <exp2> where <exp1> and <exp2> are 
strings or perl regular expressions 


<exp2> is '' by default 


Examples: 


file_rename file_new.mol2 -a _new will rename file_new.mol2 to file.mol2 


file_rename file.mol2 -a .mol2 -b .bak will rename file.mol2 to file.bak 


file_rename ssss.mol2 -a 's*' -b x -re will rename ssss.mol2 to x.mol2 


file_rename ssss.mol2 -a 's' -b x -re -s will rename ssss.mol2 to xsss.mol2 


Flags: 


-a &lt;string&gt; string to replace <br>
-b &lt;string&gt; string to insert ('' by default) <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make only a single substitution of string a (not a global one) <br>
-re&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; treat string a as a regular expression <br>
-clean&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Clean non word characters [A-Z, 0-9, - and _] from filenames <br>
-renumber&nbsp;&nbsp; Convert all numbers in filenames to have leading zeros <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; overwrite existing files <br>

Revision: 1.22.2.1.2.1 


## find_aggregate 


Identify hydrophobic aggregates of molecules (eg components of micelles) in a 
periodic system. Molecules are identified by connectivity and can contain 
multiple residues. By default, two molecules are defined as belonging to the 
same aggregate if they have carbon atoms within the distance threshold This is 
set to 4 Ang by default 


The definition of an aggregate can be changed using the -model option. The 
aggregate model can be changed to exclude polar carbons (i.e. those attached 
to O or N). This is useful for nonionic surfactants with PEG headgroups. The 
-model water option can be used to calculate water aggregates 


The -move option moves all atoms into the unit cell (molecules are often split 
when they cross the unit cell boundaries). This is useful for generating 
pictures in programs like Pymol 


A Gromacs format index file (ie atoms are indexted starting from 1) can be 
written with the -ndx flag 


Summary data is written to <first_filename>_agg.out. This file can be used 
with the unix commands 'grep AGG', 'grep COUNT' or 'grep STATS' to extract 
only information about individual aggregates, the aggregate summary or 
statistical values. 


In pdb format output, the temperature factor and SEGID columns are set to 
reflect the aggregate. In Pymol colour->spectrum->b-factor will show 
accregates 


Aggregate core atoms can be renamed to 'CX' or 'OX' using the -rename flag 


Use: find_aggregate <options> <filename> 


Flags: 


-t &lt;number&gt; Maximum distance between key atoms in aggregate (default 4 Ang) <br>
-move&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Move all atoms into the unit cell (range 0 -&gt; cell size) <br>
-i &lt;ignore_list&gt; List of molecule names to ignore. Multiple residues can be listed, separated by spaces and surrounded by brackets e.g. -i 'DAN SAM HAM'. Molecule name should be the name of the first residue in molecule. <br>
-model [carbon, nonpolar_c, water] Aggregate model. Carbon - find contacts between any carbon atoms. Nonpolar_c - find contacts betweetween carbons not bonded to N or O. Water - find water aggregates <br>
-ndx&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write Gromacs index file containing aggregates <br>
-keep&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Retain isolated molecules in calculations of average aggregate size <br>
-x &lt;number&gt; } <br>
-y &lt;number&gt; } Cell dimensions <br>
-z &lt;number&gt; } <br>
-rename&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename aggregate core atoms to 'CX' or 'OX' <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.32.2.14.2.14 


## find_amino_acids 


Identify individual amino acids within a molecule, divide into individual 
residues and provide residue names and atom names. 


Atoms are sorted into a sane order unless the -ns flag is used. 


Amino acid definitions are in SILICO_HOME/data/amino_acids_smiles.dat 


Note: Requires that hydrogen atoms are present 


Flags: 


-new&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use new routines (In development - untested) <br>
-ns&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not sort atoms <br>
-addh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force hydrogen addition <br>
-noaddh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not add hydrogens <br>

Revision: 1.8.2.1.2.18 


## find_clash 


Find clashing atoms using the find_closest_atoms routine 


CLOS segid markers are added for close atoms. These will appear in pdb format 
output file 


-clash&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use find_clash routine <br>
-ignh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore hydrogens <br>
-fix&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fix clashes by moving distance between atoms to &lt;cutoff&gt; <br>
-u&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Unbinned - do not use binned atom routines (much slower) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2013 


Revision: 1.5.2.1.2.3 


## find_close_atoms 


Flags: 


-n &lt;file&gt; Index file <br>
-r &lt;string&gt; Reference group (present in the index file) <br>
-t &lt;string&gt; Target group (optional, present in the index file) <br>
-c &lt;number&gt; Cutoff distance <br>
-i &lt;file&gt; Initial time (ps) <br>
-ts &lt;num&gt; Time step between molecules (ps) <br>
-noh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore hydrogen atoms <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use periodic boundary conditions <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print out a full report of close pairs of atoms <br>
-x &lt;number&gt; Periodic cell X axis length <br>
-y &lt;number&gt; Periodic cell Y axis length <br>
-z &lt;number&gt; Periodic cell Z axis length <br>
-O &lt;string&gt; Output file name <br>

## find_close_water 


Find water molecules that are directly hydrogen bonded to a protein 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.17.2.1.2.2 


## find_groups 


Find functional groups in a molecule using subroutine 
mol_label_functional_groups() 


Flags: 


-het&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Also print out information about aromatic rings found using mol_label_heterocycles(). <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print out internal Silico information for debugging <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


## find_ligand_contacts 


Find protein contacts for a ligand. Provide Silico atom selection strings 
suitable for use by the pacs_ligand script 


Usage example: find_ligand_contacts -n 6 -bcr 1 -nbcr 1 4buo.pdb 


This example will find protein-ligand contacts for molecule number 6 in the 
file 4buo.pdb. One backbone and one non-backbone contact per residue will be 
listed 


Flags: 


-c &lt;real&gt; Maximum contact distance (default 5.0 Ang) <br>
-n &lt;int&gt; Molecule number of ligand <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use periodic boundary conditions to measure distances (only works for cuboid cells) <br>
-all&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Find all ligand contact residues <br>
-cr &lt;int&gt; Find this number of contacts per ligand residue <br>
-pcr &lt;int&gt; Number of 'polar' contacts per ligand residue (i.e. ligand atom is O or N) <br>
-npcr &lt;int&gt; Number of 'nonpolar' ontacts per ligand residue (i.e. ligand atom is not O or N) <br>
-bcr &lt;int&gt; Number of backbone contacts per ligand residue (i.e. atoms are peptide C, N, O, CA or N) <br>
-nbcr &lt;int&gt; Number of non-backbone contacts per ligand residue (i.e. atoms are not peptide C, N, O, CA or N) <br>
-ca &lt;int&gt; Maximum number of contacts per ligand atom <br>
-min &lt;int&gt; Retain only molecules with at least this number of atoms (default 10) <br>
-split&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Format output strings by ligand residue <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC November 1999-2025 


Revision: 1.15-80-g67aad2a 


## find_rings 


Test script to find rings, planar rings and aromatic rings in a molecule. 
Temperature factor and occupancies of the output pdb file are set. The Mol2 
output file contains the rings, aromatic rings and planar rings in static 
sets. 


Flags: 


-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum depth of search for rings <br>
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out molecule structure files with rings marked <br>
-timing&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>
-debug&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; rings <br>

Revision: 1.23.2.1.2.1 


## find_similar 


Identify duplicate or similar compounds in two sets of molecules using 
Tanimoto or Euclidian comparisons 


Useage find_similar file1 file2 file3 .... 




Designed to be used to filter docking results that are ordered from best to 
worst. 


Fragments are generated using the silico fragment routines 


Tanimoto coeff = Num common fragments / Num fragments in mol1 + Num fragments 
in mol2 - Num common fragments 


Euclidian coeff = 


Flags: 


-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Scoring method (Tanimoto or Euclidian) <br>
-cut&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cutoff value for duplicate compound (default 0.80) <br>
-dup&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Find duplicates (sets cutoff to 0.999) <br>
-max&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Max fragment size <br>
-min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Min fragment size <br>
-noh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore hydrogens in fingerprints (default) <br>
-o&lt;format&gt; Output format <br>

Revision: 1.21.2.1.2.1 


## fix_cis_amides 


Crude (but effective) convert cis amides to trans. Generates structures that, 
when minimised, contain trans amides 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.4.2.1.2.3 


## formatdoc 


Print formatted comments from silico files 


Revision: 1.17.2.1.2.1 


Flags: 


-h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write HTML output to a file <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Generate subroutine descriptions <br>

## make_box 


Produce a box file '<filebase>_box.pdb'. If the input file has a unit cell 
predefined, that will be used. Otherwise, a box will be made large enough to 
enclose the first molecule. 


Unit cells are encoded in some file types, for example Gromacs, Mol2, PDB 


Flags: 


-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore existing unit cell data <br>

## make_index 


Create a an index file in Gromacs or DCD format. 


Can use a versatile atom selection language and will also write out files 
containing selected atoms 


Important: For DCD format index files see the -z option below. 


Flags: 


-z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number atoms from zero, as for DCD files (Note if you wish to use this file with catdcd, you will need to edit the resulting file to contain only a single index group and remove all lines containing square brackets <br>
-as &lt;string&gt; Use Atom Specifier (see below) <br>
-g&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create an index group containing all atoms <br>
-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create a separate index group for each atom (identified by unique atom name in each residue type) <br>
-e&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create an index group for each element <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create an index group for each residue <br>
-seg&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create an index group for each segment <br>
-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create separate index groups for water and Not_Water <br>
-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create index group for ions (residue names NA, CL, K) <br>
-wi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create index groups Water_and_Ions and Not_Water_and_Ions <br>
-set&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Create an index group for each Mol2 atom set <br>
-an &lt;atom_names&gt; Create index groups for listed atom names <br>
-rn &lt;res_names&gt; Create index groups for listed residue names <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Analyse the molecule as for a Starmaker dendrimer (create index groups for COR, GAA, GAB,...) <br>
-np&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make index of nonpolar carbon atoms (i.e. not attached to N or O). This is useful for analysing surfactants <br>
-peg&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make index of polyethylene glycol (PEG) atoms <br>
-memb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make MEMB index of membrane residues starting with the first 3 letters of standard residue names (CHL* POP* DPC*) <br>
-solu&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make SOLU index of solute residues (not water, ions or membrane residues) <br>
-cyl &lt;dist&gt; Modify MEMB residues to select only those within an XY distance &lt;dist&gt; of ligand (requires both -mem and -sol flags). For umbrella sampling simulations. <br>
-system&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make a SYSTEM index containing all atoms <br>
-write&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out a file containing atoms in index groups. The output file contains structures corresponding to each index group (except the one containing all atoms) <br>
-base&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Base name for index groups select using the atom selection language (otherwise index groups are the AS string itself)  <br>



#### Atom specifiers 


-as&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom specifier strings (multiple -as flags are allowed) <br>

Several flags are available as atom specifier shortcuts 


-back&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Backbone atoms <br>
-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA atoms <br>
-cacb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA and CB atoms <br>
-heavy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Nonhydrogen atoms <br>
-sas&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Short specifier [ca, cacb, back, back_disulfide, heavy] <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Read atom specifiers from a file <br>

#### Atom specifier examples: 


ANAME:CA&nbsp;&nbsp;&nbsp;&nbsp; All atoms called 'CA'. <br>
ANAME:CA,CB,CG Returns all atoms called 'CA' or 'CB' or 'CD'. <br>
ANAME:CA,RESNAME:TRP All atoms called 'CA' in all residues called TRP <br>
ANAME:CA,SUBID:4 All atoms called 'CA' in residue number 4 <br>
ELEMENT:!H All nonhydrogen atoms <br>
SEGID:PROT,SEGID:LIG All atoms with the SEGID set to PROT or LIG <br>

Successive atom specifications can be made. Each is separated by a '|' 


ANAME:CB|ANAME:CA|ANAME:CD Returns all atoms called 'CB', 'CA', 'CD'. <br>
ANAME:CA,RESNUM:1|ANAME:CA,RESNUM:4 Returns 'CA' atoms from residue 1 and 4. <br>

Atom specifiers are case sensitive. 


Revision: 1.30.2.1.2.33 


## measure_dist 


Find the shortest distance between a specified atom and a specified residue 


Useful for investigating GPCR-ligand interactions 


Use options -i and -l <ligand_name> for output containing both ligand and 
protein in the same structure (e.g. Glide induced fit output) 


Flags: 


-p &lt;protfile&gt; Protein file name (if using separate file) <br>
-r &lt;resnum&gt; Protein residue number <br>
-la &lt;name&gt; Ligand atom name <br>
-lr &lt;resnum&gt; Ligand residue name <br>
-ln &lt;num&gt; Ligand atom number <br>
-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Input structure contains both protein and ligand (e.g. From Glide IFD docking) <br>
-f &lt;dist&gt; Filter out ligands with n_dist greater than this value. 0 = ignore. <br>
-s &lt;string&gt; Skip ligands with names starting with 'string' <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2 


Created: DKC 2010 and later 


## Micelle 


Build a spherical micelle using based on regular distributions of points on a 
sphere 


For a similar program with random placement of molecules see vesicle_builder 


Use: micelle -n <num_lipids> <lipid> 


Expects the <lipid> to be aligned along the +Z axis with the hydrophobic end 
located at the origin (0,0,0). Use mol_rotrans with the -alignz option to do 
this. The -offset option can be used to prevent or reduce clashes caused by 
overlapping atoms at the origin. 


Not all numbers <num> are possible. Use the -l flag to list available numbers 


Flags: 


-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List available packing numbers and exit <br>
-n &lt;num&gt; Number of lipids to pack <br>
-offset &lt;dist&gt; Distance to offset lipid along z axis <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.5 


Created: DKC March 2017 and later. 


## mmod_csearch 


Run a conformational search using Macromodel 


Flags: 


-ff &lt;value&gt; Forcefield [opls2005, opls2001, amber94, amber, mm3, mm2] (default opls3) <br>
-s &lt;solvent&gt; GBSA solvent [water chloroform octanol none] (default water) <br>
-cw&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Run search in chloroform and water solvents <br>
-method &lt;value&gt; Conformational search method [mcmm, spmc, lmc2, lmcs, lmcs-mix, lmcs-shake, lmcs-mode (Default mcmm) <br>
-flap&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use FLAP option for ring searchers in mcmm and low-mode searches (not implemented yet) <br>
-ns &lt;steps&gt; Number of conformational search steps (default 5000) <br>
-min &lt;steps&gt; Maximum number of steps for minimisation <br>
-grad &lt;grad&gt; Gradient cutoff for minimisation <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Overwrite existing output files <br>
-maxtok &lt;num&gt; Maximum number of Schrodinger tokens to use <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.8 


## mmod_min 


Minimise a molecule using Macromodel 


Flags: 


-ff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Forcefield [opls2005, opls2001, amber94, amber, mm3, mm2] (default opls3) <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; GBSA solvent [water chloroform octanol none] (default water) <br>
-min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of steps for minimisation <br>
-grad&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gradient cutoff for minimisation <br>
-minh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Minimise hydrogens only <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.6.6.5 


## schrod_qikprop 


Run Schrodinger Qikprop 


Flags: 


-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Overwrite existing output files <br>
-maxtok &lt;num&gt; Maximum number of Schrodinger tokens to use <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.1 


## mol2seq 


Extract amino acid sequence from molecule 


Output file <filebase>.xxx 


Flags 


-ih&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include HETATMs <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Combine all sequences into a single file (all.seq) <br>
-o &lt;format&gt; Output format <br>

Revision: 1.27.2.1.2.1 


## mol2split 


Fast split program to divide a mol2 file into smaller files without parsing 
the file 


Default is 100 structures per file 


Flags: 


-s &lt;style&gt; Output style [numbered, fnumbered] (default: fnumbered) <br>
-n &lt;number&gt; Number of structures in each output file (default: 100) <br>
-d &lt;dir&gt; Output directory (default: working directory) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.17.2.1.2.2 


## mol_add_h 


Add hydrogens to a molecule 


By default, a protonation state is produced which approximates physiological 
state. Using the -v flag will fill all valences. ie Will add one hydrogen to a 
carboxylic acid and 3 hydrogens to ammonia 


Adds only polar hydrogens if the 'polar' flag is used. Adds hydrogens to 
carbon only if 'nonpolar' flag is used 


Flags: 


-polar&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add only polar hydrogens <br>
-nonpolar&nbsp;&nbsp; Add only nonpolar hydrogens (hydrogens on carbon) <br>
-v&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fill valence <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; File containing atom names and connectivities ($SILICO_HOME/data/amino_acid_atoms.dat) <br>
-ns&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not sort atoms <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print atoms with formal charges <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001-3 


Revision: 1.32.2.1.2.4 


## mol_add_lp 


Add lone pairs to a molecule 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2007 


Revision: 1.7.2.1.2.2 


## mol_add_segid 


Add a SEGID and write a pdb file 


Flags: 


-segid &lt;SEGID&gt; <br>
-i &lt;input_format&gt; Input format, any or pdb (default any) <br>
-o &lt;output_format&gt; Output format, all silico formats (default pdb) <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delete disordered (ALT) atoms <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; use SINGLE option <br>
-prefix&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add 'new_' to start of output file and do not add '_new' to end <br>
-debug&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.4.2.1.2.1 


## mol_amides 


Constructs a plot of dihedral angles w and w' for secondary amides 


Flags: 


-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print a hardcopy instead of producing an output file <br>
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;format&gt; Output file format (default: ps) <br>
-g&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to grace executable <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: BPR September - October 2004 


Revision: 1.25.2.1.2.3 


## mol_backbone 


Delete all non-backbone atoms from a molecule. 


Uses atom labels to determine which atoms should be kept. Default set is {N, 
CA, C, O}. This can be further reduced to CA atoms only by use of the -ca 
flag, as below. 


Flags: 


-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Keep only CA atoms <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.16.2.1.2.3 


## mol_boltz_pop 


Calculate the relative Boltzmann populations for an ensemble of molecules. 
Assumes energies are in kJ/mol and that the first structure is the lowest 
energy. 


Flags: 


-temp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Temperature (default 310K, 37 C) <br>
-ef Energy field (default 'Relative_Potential_Energy-OPLS-2005') <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2 


## mol_bond 


Test bond creation 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC November 1999-2025 


Revision: 1.1.2.2.2.9 


## mol_centre 


Translate a molecule centre of cell (default) or origin 


Flags: 


-origin&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate centre to origin <br>
-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Centre based on each molecule's largest fragment <br>

Geometry Flags 


-x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>
-y&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } Dimensions of the box (Angstroms) <br>
-z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 


Revision: 1.20.2.1.2.6 


## mol_centroid 


Calculate the centroid of a molecule and add a dummy atom at that point 


Flags: 


-x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>
-y&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } Dimensions of the box (Angstroms) <br>
-z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 


Revision: 1.15-61-gef33873 


## mol_characterise 


Calculate molecular weight, molecule extents and other molecular properties 


Properties 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of atoms  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of bonds  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of rotatable bonds (see subroutine mol_count_rot_bonds)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Molecular weight  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of C, H, N, O, P, S, halogen and other atoms.  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of rings up to size 10 (Note that the number of rings is not quite the way a chemist would see it - eg Naphthalene has 3 rings)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of planar rings  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of H-bond donors (see subroutine mol_find_donors_acceptors)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of H-bond acceptors  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Molecule name  <br>

Adds hydrogens first by default 


Data is written to a file '<filebase>_mc.out' in tab delimited format and into 
SDF data fields if -sdf flag is set 


Flags 


-sdf&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out an sdf file containing molecule structure and calculated data (on by default) <br>
-hadd&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add hydrogens (on by default) <br>
-t&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write data to tab delimited text file <br>
-replace&nbsp;&nbsp;&nbsp;&nbsp; Overwrite original sdf file (only if -sdf is set) <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Overwrite preexisting output files <br>

Revision: 1.20.2.1.2.2 


## mol_charge 


Find the total partial and formal charges on molecle 


Flags: 


-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Calculate formal charges <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print atoms with formal charges <br>
-df&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Die if formal charge is not zero <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Provide a charge breakdown by residue <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Provide a charge breakdown by segment <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.15-59-g89487b8 


## mol_check 


Run a series of sanity checks on a molecule 


Calls ensemble_check (this principally checks the integrity of the internal 
Silico data structures), mol_check_atom_overlap (to find badly overlapping 
atoms), mol_check_valences (to find atoms with an incorrect number of 
connected atoms) 


Writes out a mol2 file with subsets ERROR, OVERLAP, BONDLENGTH and VALENCE 
containing any atoms with errors 


It would be desirable to check for poor bond lengths and angles as well 


Flags 


-ce&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check atom elements (default on) <br>
-co&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check for atom overlap (default off) <br>
-cv&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check atom valences (default on) <br>
-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check for aromatic bonds in non-aromatic systems (default on) <br>
-cb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check bondlengths (default on) <br>
-cr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check for distorted aromatic rings (default on) <br>
-amide&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check for cis-amides and distorted trans amides (default on) <br>
-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Run all checks <br>
-noconnect Do not create connection table <br>
-nofile&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not write an output .mol2 file for each input file <br>

Created: DKC 2000 


Revision: 1.15-79-g6958a7c 


## mol_combine 


Combine separate molecules into a single molecule 


By default combines all the molecules on the command line in to one single 
molecule. Each molecule is given a separate SEGID 


If the -p option is used to specify a 'parent' molecule, then each molecule 
specified on the command line will be combined separately with the parent 
molecule. This is useful for combining many ligands with a single parent 
protein to produce complexes. 


The -glide option will combine glide _pv.maegz files to produce multiple 
receptor/ligand structures in individual files. If the -o pdb option is chosen 
then the ligand will have HETATM records to be compatible with ligplot. 


Should work with gromacs .itp files (not well tested yet...) 


See also: mol_unsplit - combines molecules into a multi-molecule file 




Flags: 


-o &lt;format&gt; Output file format <br>
-p &lt;parent protein&gt; Combine all molecules with parent (usually protein) molecule <br>
-glide&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Combine all molecules with first molecule in file. Use with glide .pv or .raw files <br>
-ra&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Renumber atoms <br>
-dr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Don't renumber residues i <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 and onwards 


Revision: 1.29.2.4.2.9 


## mol_del_atoms 


Delete atoms from a file 


Atoms are selected on the basis of their names, their elements or their 
residue names. 


Name, element and residue name options accept space separated values, which 
are ORed. For example, -e C,N would mean "delete carbons or nitrogens". 


The * wildcard can be used in atom or residue names, but it must be enclosed 
in quotes to escape the shell. For example, mol_del_atoms -a N'*' file.mol2 
will delete atoms whose names start with N in all molecules in file.mol2. 


The various criteria are ANDed. For example, -r HOH -e O,H would mean "delete 
all atoms whose residue name is HOH and which are oxygens or hydrogens". 


If any given criterion is left blank, any value for that criterion is 
considered acceptable. 


Use silico atom specifiers given in on the command line using the -as flag or 
in a file using the -asf flag 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom names (space separated list) <br>
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom numbers (space separated list, atoms are numbered starting at 1) <br>
-e&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom elements (space separated list) <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Residue names (space separated list) <br>
-as &lt;string&gt; Use Atom Specifier <br>
-asf &lt;filename&gt; Filename containing atom specifier <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Combine all input files to a single output file <br>
-q&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Speed up deletion if all molecules in file are identical  <br>



#### Atom specifiers 


-as&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom specifier strings (multiple -as flags are allowed) <br>

Several flags are available as atom specifier shortcuts 


-back&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Backbone atoms <br>
-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA atoms <br>
-cacb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA and CB atoms <br>
-heavy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Nonhydrogen atoms <br>
-sas&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Short specifier [ca, cacb, back, back_disulfide, heavy] <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Read atom specifiers from a file <br>

#### Atom specifier examples: 


ANAME:CA&nbsp;&nbsp;&nbsp;&nbsp; All atoms called 'CA'. <br>
ANAME:CA,CB,CG Returns all atoms called 'CA' or 'CB' or 'CD'. <br>
ANAME:CA,RESNAME:TRP All atoms called 'CA' in all residues called TRP <br>
ANAME:CA,SUBID:4 All atoms called 'CA' in residue number 4 <br>
ELEMENT:!H All nonhydrogen atoms <br>
SEGID:PROT,SEGID:LIG All atoms with the SEGID set to PROT or LIG <br>

Successive atom specifications can be made. Each is separated by a '|' 


ANAME:CB|ANAME:CA|ANAME:CD Returns all atoms called 'CB', 'CA', 'CD'. <br>
ANAME:CA,RESNUM:1|ANAME:CA,RESNUM:4 Returns 'CA' atoms from residue 1 and 4. <br>

Atom specifiers are case sensitive. 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.21.2.1.2.8 


## mol_del_dummy 


Delete all dummy atoms from a file (eg lone pairs). 


Output file <filebase>_nodu.<ext> 


Flags: 


-con&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force regeneration of connection table <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2/1/2000 


Revision: 1.18.2.1.2.1 


## mol_del_h 


Delete all hydrogens from a file. 


Output file <filebase>_noh.<ext> 


Flags: 


-res&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delete hydrogens on a particular residue name <br>
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delete all nonpolar hydrogens <br>
-nsp2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delete all nonpolar hydrogens except sp2 hydrogens (aromatics, aldehydes, etc) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2/1/2000 


Revision: 1.25.2.1.2.3 


## mol_del_solvent 


Delete solvent molecules from a file (or alternatively delete nonsolvent 
molecules) 


Writes out _dry file containing unsolvated molecules or _sol file containing 
solvent 


Flags: 


-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Negate. Write out solvent instead of nonsolvent <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force deletion of water and other solvent molecules <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Combine all input files to a single output file <br>
-q&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Speed up deletion if all molecules in file are identical <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.7 


## mol_divide 


Separate a multi-molecule file (eg Tripos mol2 or Schrodinger mae) containing 
into individual files, each containing one molecule. 


Files are put into a specified directory. This is named <filebase>.dir if a 
name is not provided. Each molecule is renamed with an Insight-safe name. The 
default behaviour renames the output file to match the molecule name. Other 
file naming styles can be selected using the -s flag 


To separate a single molecule into multiple structures see 'mol_split' 


For really big files (multiple thousands of structures) consider using 
sdfsplit, mol2split or pdbsplit which do not parse the file and are much 
faster 


For convenience the script makes a pymol load script 'load.pml' in the output 
directory. Run this script in pymol to load all the files 


For more control over splitting PDB files by chain, waters etc see 'pdbsplit' 


Flags: 


-style&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Output style for file names. molname: molecule name (not checked for duplicates!). molname_i 'insight_safe' molecule name. numbered: numbered, no leading zeros, fnumbered: numbered, leading zeros <br>
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of structures per file (default 1) <br>
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Output format <br>
-force&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Overwrite existing output <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make pymol load script 'load.pml' in output directory <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Directory name base <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC October 2000 and on... 


Revision: 1.44.2.1.2.9 


## mol_ensemble_average 


Given a number of molecules of the same composition as input, write out a 
molecule where the position of each atom is averaged over the whole ensemble. 


Flags: 


-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Consider the largest fragment only <br>
-t&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate the molecule (or largest fragment if -l) to centre of mass <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: BPR 2007 


Revision: 1.7.2.1.2.1 


## mol_extents 


Calculate the maximum and minimum X, Y and Z coordinates of a molecule and the 
centre point. 


Values are written to STDOUT 


Flags: 


none&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <br>

Revision: 1.18.2.1.2.1 


## mol_filter 


Filter a set of molecules by SDF property and/or molecular weight. 


The 'noparse' option can greatly speed up proscessing of SDF and some other 
formats by avoiding processing of most atom and bond information. This is not 
compatible with the -fg option. 


Can be used to retain molecules with a unique SDF_CODE 


Flags: 


-p &lt;name&gt; SDF property name <br>
-max &lt;val&gt; Property maximum value <br>
-min &lt;val&gt; Property minimum value <br>
-mwmax &lt;val&gt; MW maximum value <br>
-mwmin &lt;val&gt; MW minimum value <br>
-hamax &lt;val&gt; Heavy atom maximum value <br>
-hamin &lt;val&gt; Heavy atom minimum value <br>
-druglike&nbsp;&nbsp; Select druglike compounds (max MW 500, min MW 150) <br>
-leadlike&nbsp;&nbsp; Select leadlike compounds (max MW 400, min MW 75) <br>
-fragmentlike Select fragmentlike compounds (MW 300, min MW 60) <br>
-u&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Retain only a single molecule foreach unique value of this SDF field. e.g. retain only a single molecule for each molecule NAME <br>
-fg &lt;list_of_groups,_space_delimited&gt; Retain only compounds containing these functional groups. Functional group names are determined by mol_label_functional_groups <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;string&gt; Select compounds containing &lt;string&gt; in the id field (case insensitive). The id field is set using the -id flag <br>
-cw&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The -c flag must only match a whole word <br>
-list &lt;filename&gt; Select compounds matching a list provided in a file. File must contain one string to be matched per line <br>
-id &lt;id&gt; SDF field for list filter <br>
-v&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Verbose. Print information about discarded compounds <br>
-np&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use noparse option (much faster, but currently only available for sdf files) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2003 


Revision: 1.32.2.3.2.5 


## mol_flatten 


Squash molecules flat into the XY plane 


Revision: 1.18.2.1.2.1 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


## mol_hydrogen_bonds 


Find all hydrogen bonds in a system. 


Uses the hydrogen bond definition of McDonald and Thornton (J. Mol. Biol. 
1994, 238, 777-793) which specifies maximum distances and angles for A...H-D 
and A..D. Note that to increase the H-bond distance, you must increase both 
the A..H-D and A..D distances. 


Output molecule file contains structures with hydrogen bonds are denoted using 
0-order bonds. These can be visualised in Maestro or Pymol (etc). 


The default output filename is derived from the first file. This can be 
changed using the -O flag. 




Output file lists the following information 


Energy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Energy value derived from input file <br>
Backbone-backbone H-bonds (between NH and C=O) <br>
3-10&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Peptide backbone H-bonds in 3-10 helix <br>
alpha&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Peptide backbone H-bonds in alpha helix <br>
NumDon&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total number of H-bond donors <br>
NumAcc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Total number of H-bond acceptors <br>
NDon&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of nitrogen H-bond donors <br>
ODon&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of oxygen H-bond donors <br>
NAcc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of nitrogen H-bond acceptors <br>
OAcc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of oxygen H-bond acceptors <br>
UnsatDon&nbsp;&nbsp;&nbsp;&nbsp; Number of unsatisfied H-bond donors (zero hydrogen bonds to this atom) <br>
UnsatDon&nbsp;&nbsp;&nbsp;&nbsp; Number of unsatisfied H-bond acceptors (zero hydrogen bonds to this atom) <br>
NumChgs&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of formally charged atoms (carboxylate counts as 1,etc) <br>
ChgReinf1&nbsp;&nbsp; Number of charge reinforced H-bonds with either donor or acceptor charged (delocalisation on carboxylate groups and guanidines is considered). <br>
ChgReinf2&nbsp;&nbsp; Number of charge reinforced H-bonds with both donor and acceptor charged (delocalisation on carboxylate groups and guanidines is considered). <br>

Flags: 


Hydrogen bond parameters 


-d &lt;val&gt; Maximum Donor-Acceptor distance (default 3.9 Ang) <br>
-h &lt;val&gt; Maximum Hydrogen-Acceptor distance (default 2.5 Ang) <br>
-a &lt;val&gt; Minimum Donor-Hydrogen-Acceptor angle (default 90 deg) <br>
-b &lt;val&gt; Minimum Hydrogen-Acceptor-Substituent angle (default 90 deg) <br>

Timestep parameters 


-ts &lt;val&gt; Timestep between files (ps) <br>
-i &lt;val&gt; Time at first file (ps) <br>

Atom/molecule selection options 


-wat&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include water molecules <br>

Input file options 


-copy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Copy data from first molecule to subsequent molecule. This is good for MD trajectories <br>

Output file options 


-nhb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not create 0-order bonds for all hydrogen bonds and series of PDB files <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: BPR 2005 


Revision: 1.30.2.1.2.15 


## mol_label_fg 


Test script to test the subroutine 'mol_label_functional_group' 


Revision: 1.15-89-gc2af617 


Flags: 


-aa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Label aa backbone <br>
-het&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Label heterocycles <br>
-seg&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Set SEGID to functional group name for debugging <br>

## mol_mw 


Calculate molecular weight and molecule extents 


MW data is calculated for the parent molecule (molecule with most atoms in 
file) 


Values are printed to standard output and added to output file as SDF_DATA 


Revision: 1.26.2.2.2.3 


Flags: 


-addh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add hydrogens <br>
-v&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fill valence <br>
-print&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print data to a file 'mw.txt' <br>
-min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Minimum molecular weight to output (ie skip molecules below this mass) <br>
-max&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum molecular weight to output (ie skip molecules above this mass) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


## mol_orient 


Orient molecules so that their axis of polarity is in the XZ (or XY) plane 
Useful for building bilayers and other structures 


By deffault molecules are aligned along the Y axis with the nonpolar centre at 
the origin and with the polar centre in the +Y direction. This is sutable for 
bulding a bilayer in the XZ plane (the Silico default). 


Use the -xy flag to orient the molecule along the Z axis for building a 
bilayer in the XY plane. 


Use: mol_orient <molecules...> 


Flags: 


-xy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Build bilayer in the xy plane <br>
-add&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add dummy atoms at polar and nonpolar centres <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.15-7-gf32c0c6 


Created: DKC November 1999 and later. 


## mol_printout 


Print out internal Silico data for a molecule 


Flags: 


-all&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print all atom recors (otherwise only 10 are printed) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 and later 


Revision: 1.1.2.1 


## mol_rama 


Constructs a Ramachandran plot for a molecule containing Alpha Amino acids 


Plotting is performed by 'xmgrace' 


Flags: 


-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Output format (default: PostScript) <br>
-print&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print a hardcopy (default: off) <br>
-label&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Label points with residue name <br>
-residue&nbsp;&nbsp;&nbsp;&nbsp; Make one plot for each residue <br>
-mol&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make one plot for each molecule in file. Use if file contains different molecular structures <br>
-ramol&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make a separate .ramol file with all dihedral angles for one molecule on a single line <br>
-energy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include energy in ramol output file <br>

Created: BPR September - October 2004 


Revision: 1.15-60-g738dd58 


## mol_rename 


Rename molecules 


Can also: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;change the molecule name to the filename  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;change the molecule filename to the molecule name  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rename the molecule using a specified SDF data field  <br>

Using both the -r<sdf_datafield> and -c flags can be used to change filename 
to the specified SDF data field. 


Can generate 'insight_safe' names. i.e. So that they do not contain spaces, 
punctuation, start with an underscore or a digit and are of limited length. 


Flags: 


### Options that modify the molecule name 


-b &lt;name&gt; Set the molecule name to &lt;name&gt; (can be also used as a molecule base name) <br>
-a &lt;suffix&gt; Add &lt;suffix&gt; to the end of the existing molecule name <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Change molecule name to filename <br>
-r &lt;field_name&gt; Rename molecules using the &lt;field_name&gt; SDF Data field. <br>
-smiles&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Renamie molecules using SMILES string  <br>

These flags modify the molecule name in the CT record of SDF files 


These three flags can be combined. The changes will be concatenated in the 
order -b, -f, -r 


Options the modify the existing molecule name 


-safe&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use insight-safe names (ie that do not contain spaces, punctuation, start with an underscore or a digit and are of limited length) <br>
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add a number to the end of the name <br>

These flags can be specified in conjunction with the flags above 


### Options that modify molecule data (SDF data) 


-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Transfer the molecule name to the SDF Data fields 'NAME' and 'title' <br>
-sdfield &lt;field_name&gt; Change specified SDF field to molecule name <br>
-sddel &lt;"field_field"&gt; Delete SDF field  <br>

### Options that make multiple changes to the molecule name 


-g &lt;base&gt; Generate new name using &lt;base&gt;, add a number and change SDF Name field (same as -b, -n, -s)  <br>

### Options that modify the file name 


-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Change output filename to name of first molecule  <br>

### Other options 


-k&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Keep the same filename, overwriting the input file (possibly dangerous) <br>
-np&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use noparse option (much faster, but currently only available for sdf files)  <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000-2014 


Revision: 1.37.2.5.2.9 


## mol_renumber 


Renumber and/or rename atoms and/or residues in a molecule 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;Starting atom number&gt; Renumber atoms starting from this number <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;Starting residue number&gt; Renumber residues starting from this number <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;Starting chain letter&gt; Relabel chains starting from this letter atom they are connected to <br>
-rename&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename atoms <br>
-method &lt;val&gt; Atom renaming method (consecutive, element, element-h). Consecutive: atom numbers are increased by 1 for each atom. Element: Atom numbers increase per element. Element-H: Atom numbers increase per element with hydrogens named after the parent atom <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.26.2.3.2.3 


## mol_rescale_bonds 


Rescale a molecule so that the carbon-carbon bonds have a reasonable average 
bond length. 


This defaults to 1.5 Angstroms, however, it can be adjusted through use of the 
-l flag. Alternatively, a scaling factor can be used by means of the -f flag. 


This is useful to clean up files that have come out of Isis databases before 
they are minimised by some other program (eg Insight). 


Flags: 


-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;factor&gt; Scaling factor (not to be used with -l) <br>
-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;length&gt; Target C-C bond length (not to be used with -f) <br>
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;format&gt; Output format <br>
-split&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split multiple molecule files into a separate directory. Each molecule is renamed with an "Insight safe" name <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC Thu Sep 21 13:19:53 EST 2000 


Revision: 1.22.2.1.2.4 


## mol_rot 


rotate a molecule about any vector 


F 


-x &lt;number&gt; } <br>
-y &lt;number&gt; } Vector to rotate around (assumed to pass thru origin) <br>
-z &lt;number&gt; } <br>
-test&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Testing routine: rotate the first molecule only about the axis by 60 degrees <br>
-a &lt;number&gt; Angle to rotate molecule through <br>
-random &lt;number&gt; Generate this number of randomly rotated molecules <br>
-maximise&nbsp;&nbsp; Approximately maximise the extents of the molecule along the X and Y axes (ie align the major axis along the X axis and the medium axis along the Y axis). <br>

SF 


Revision: 1.25.2.1.2.1 


## mol_rot_bond 


Set the torsion angle between four atoms to a specified value. 


The four atoms are not necessarily bonded to each other, however it makes more 
chemical sense if they are (provided the middle two are not in a ring). If the 
middle two atoms are in a ring, nothing will be done. 


Using the -inc and -steps flag will generate additional strucutres incremented 
by angle 'inc' 


Flags: 


-a &lt;number&gt; Atom A number <br>
-b &lt;number&gt; Atom B number <br>
-c &lt;number&gt; Atom C number <br>
-d &lt;number&gt; Atom D number <br>
-ang &lt;number&gt; Desired (or initial) torsion angle (degrees) <br>
-inc &lt;number&gt; Angle increment (optional) <br>
-s &lt;number&gt; Number of steps to increment (optional) <br>

Revision: 1.16.2.1.2.5 


## mol_rotrans 


Apply rotations and/or translations to a file 


Translation is applied first, followed by rotations about the X, Y and/or Z 
(in that order) 


Flags: 


-x &lt;angle&gt; } <br>
-y &lt;angle&gt; } X, Y and Z rotation angles in degrees <br>
-z &lt;angle&gt; } <br>
-a &lt;dist&gt; } <br>
-b &lt;dist&gt; } X, Y and Z translation distances in Angstroms <br>
-c &lt;dist&gt; } <br>
-alignz&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate specified atom (default atom 1) of molecule to the origin and align the molecule along the Z axis. No other rotations or transliations are performed. <br>
-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;atomnum&gt; Atom to be moved to the origin with alignz flag (default atom 1) <br>
-invert&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make mirror image by flipping the sign of Z coordinate <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.17.2.1.2.4 


## mol_scale 


Scale a molecule by a factor 


Flags: 


-f &lt;number&gt; Scale factor <br>
-o &lt;format&gt; Output file format (default: input format) <br>
-nocell&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not scale unit cell <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.14.2.1.2.2 


## mol_segment 


Split all molecules in file to separate molecules (based on connectivities) 
and recombine them in to a single molecule. Each molecule is placed in a 
separate segment (M001 ... MXXX). 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 


Revision: 1.19.2.1.2.1 


## Mol_set_chain 


Set the PDB chain of a molecule 


Flags: 


-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Chain (default A) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 


Revision: 1.1.2.1 


## mol_smiles 


Generate a SMILES string for a molecule 


Flags: 


-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Generate smiles of parent molecule (ie remove salts, etc) <br>
-b&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include explicit bond orders <br>
-h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include explicit hydrogen atoms <br>
-k&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use Kekule bonds and non-aromatic atom symbols <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use ChemAxon molconvert to generate canonical Smiles rather than Silico <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.9.2.3.2.1 


## mol_solvate 


Make a solvated box around a molecule with a density 




Only one file may be supplied as an argument. 


Default water residue name is HOH with atoms labelled OH2, H1, H2. Using the 
-amber flag produces residue name WAT with atom names O, H1, H2 Using the 
-gromacs flag produces residue name SOL with atom names OW, HW1, HW2 


Usage: mol_solvate <file> [<flags>] 


Flags: 


Solute flags 


-t&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate solute molecule to coordinate origin <br>

Solvent flags 


-n &lt;number&gt; Add this number of solvent molecules <br>
-s &lt;file&gt; Solvent molecule file <br>
-d &lt;number&gt; Required density (default 1 g/mL) <br>
-r &lt;resname&gt; Specified name for water molecules (overrides standard water names) <br>
-amber&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Give water molecules AMBER names (resname WAT, Oxygen O, hydrogens H1, H2), and also adds TER after each water molecule <br>
-gromacs&nbsp;&nbsp;&nbsp;&nbsp; Give water molecules GROMACS names (resname SOL, Oxygen OW hydrogens HW1, HW2), <br>

Geometry flags 


-x &lt;number&gt; } <br>
-y &lt;number&gt; } Dimensions of the box in each direction <br>
-z &lt;number&gt; } <br>
-bilayer &lt;number&gt; Solvate bilayer. ie Don't put solvent molecules within &lt;number&gt; A of bilayer plane. Default is XZ. XY plane if -xy flag is used.  <br>
-xy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Solvate bilayer in the xy plane <br>
-plane &lt;thick&gt; Construct a plane of molecules in the centre of the cell of thicknes &lt;thick&gt;. By default the XZ plane is constructed unless the -xy flag is used. <br>
-sphere &lt;radius&gt; Add solvent molecules in sphere of radius &lt;radius&gt; <br>
-invert&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Invert the sense of all 3d structures <br>

Packing flags 


-si &lt;factor&gt; Initial scale factor for atom-atom packing distance (default 0.8) <br>
-sm &lt;factor&gt; Minimum scale factor for atom-atom packing distance (default 0.4) <br>
-incr &lt;number&gt; Number of unsuccesful molecule packing trials before reducing the scale factor. Default 100 <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.50.2.1.2.9 


## mol_sort 


Script to reorder the atoms in a file 


Note. -sa and -pa flags can be combined usefully 


Created: DKC 2000 


Flags: 


-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rearrange residue order (takes list of residue numbers from command line) <br>
-sa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms alphabetically (ignoring all other fields) <br>
-pa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms using 'pdb_molecule_sort_atoms' (considers residue, chain and SEGID) <br>
-sr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort by residue name <br>
-cm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Canonical SMILES sort of all atoms in molecule <br>
-rename&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename atoms after sorting <br>
-fix&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fix molecule chains and SEGID records <br>
-fixs&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fix substructure IDs (make each chain a separate substructure) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.22.2.5.2.2 


## mol_split 


Split a molecule file into separate molecules based on connectivities. 


Molecules are named using the RESIDUE name of the first residue by default 


To separate a multi-molecule file into individual structures see 'mol_divide' 


Known Bugs -d option does not work properly with pdb files 


Flags: 


-r &lt;RESNAME&gt; Retain only molecules with this residue name <br>
-min &lt;NUM&gt; Retain only molecules with at least this number of atoms <br>
-max &lt;NUM&gt; Retain only molecules with this number of atoms or fewer <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Name molecules by by SEGID <br>
-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Keep each molecule's largest fragment only <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write each output molecule as a separate file <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Renumber residues of split molecules starting from 1 <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 and on (and on) 


Revision: 1.36.2.1.2.6 


## mol_split_segid 


Split a molecule file into separate molecules based on SEGID. Molecules are 
written out to a single file. Molecules with no defined SEGID are assigned to 
the SEGID 'NONE'. 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2002 


Revision: 1.22.2.1.2.2 


## mol_transfer_coords 


Transfers X, Y and Z coordinates from file1 to file2 


Can be useful to apply crystal complex positions to a correctly 


parameterised file 


Use: mol_transfer_coords <file1> <file2> 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2.2.1 


Created: DKC November 1999 and later. 


## mol_unsplit 


Combine molecule files into a multiple molecule file 


To combine molecules into a single molecule see 'mol_combine' 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 and on (and on) 


Revision: 1.1.2.2 


## mol_wrap_cell 


Wrap all molecules back in to a unit cell. 


By default the unit cell starts at 0, 0, 0 and the largest fragment is moved 
to the centre of this cell 


The system around a given atom using the -a flag 


Useful for molecular dynamics output where some molecules have wandered out of 
the unit cell 


Flags: 


-ignore&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ignore unit cell dimensions in individual files <br>
-x &lt;number&gt; Default unit cell's X dimension <br>
-y &lt;number&gt; Default unit cell's Y dimension <br>
-z &lt;number&gt; Default unit cell's Z dimension <br>
-i &lt;file&gt; Index file <br>
-a &lt;number&gt; Atom number to centre the system around <br>
-g &lt;string&gt; Name of index group to centre the system around <br>
-t&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate the centre of largest molecule to (0,0,0) <br>
-l&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate bottom, left, (0,0,0) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 


Revision: 1.31.2.1.2.4 


## n_dist 


Find the shortest distance between a charged nitrogen in the ligand and 
carboxylate atoms of a specified residue (residue must contain a carboxylic 
acid - e.g. Asp, Glu). 


Useful for investigating GPCR-ligand interactions 


Use options -i and -l <ligand_name> for output containing both ligand and 
protein in the same structure (e.g. Glide induced fit output) 


Flags: 


-r &lt;resnum&gt; Residue number of acidic residue in protein <br>
-p &lt;protfile&gt; Protein file name <br>
-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Input structure contains both protein and ligand (e.g. From Glide IFD docking) <br>
-f &lt;dist&gt; Filter out ligands with n_dist greater than this value. 0 = ignore. Ligands without a charged nitrogen are also skipped. <br>
-s &lt;string&gt; Skip ligands with names starting with 'string' <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.6.2.3.2.3 


Created: DKC 2022 


## pdb_chain_split 


Split a pdb file into chains keeping any ligands with that chain 


Use: pdb_chain_split <file1> <file2> ... 


Flags: 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2 


Created: DKC November 1999 and later. 


## pdb_rename_hydrogens 


Rename hydrogen atoms so that they have the correct PDB nomenclature. Using 
the -charmm flag will produce charmm27 atom names. Using the -cyana flag will 
produce cyana2 names. 


Note: Particular attention must be paid to the delta-carbon of isoleucine 
residues, which is also renamed. The -charmm flag will name ILE CD as CD. The 
-cyana flag will name ILE CD as CD1. 


C-terminus and N-terminus hydrogen names are currently not generated. 


Flags: 


-charmm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Generate charmm27 atom names <br>
--charmm-atom-names <br>
-cyana&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Generate cyana2 atom names <br>
--cyana-atom-names <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;filename&gt; Filename containing atom names and connectivities <br>
--datafile=&lt;filename&gt; <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delete any hydrogens that can not be given a name <br>
--delete-unknown-h <br>
-o&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;format&gt; Output file format (default PDB) <br>
--output-file-format=&lt;format&gt; <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


David K. Chalmers, 3 February 2000 


Revision: 1.17.2.1.2.3 


## pdbsplit 


Split a PDB file into smaller pieces based on TER or MODEL records without 
parsing the file 


Use: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pdbsplit file.pdb -s fnumbered produces file_00001.pdb, file_00002.pdb, etc  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;pdbsplit file.pdb -s numbered produces file_1.pdb, file_2.pdb  <br>

Flags: 


-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Filename output style <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;dirname&gt; Output directory <br>
-end&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split on END records <br>
-endmdl&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split on ENDMDL records <br>
-ter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split on TER records <br>
-chain&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split at end of each chain/SEGID. Name output files by chain <br>
-all&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split at all of the above (default but unset by selecting one of the above) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.20.2.1.2.1 


## plot_dihedrals 


Wrapper script to plot dihedral potentials using gnuplot. 


Requires:&nbsp;&nbsp; Either an input file, or a number of dihedrals to plot <br>
Returns:&nbsp;&nbsp;&nbsp;&nbsp; Plot <br>

Revision: 1.20.2.1.2.1 


## radius_of_gyration 


Calculates the radius of gyration for a series of molecules 


Flags: 


-ts &lt;number&gt; Timestep between files (ps) <br>
-i &lt;number&gt; Time at initial file (ps) <br>
-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use all fragments (not just the largest one) to calculate radius of gyration <br>
-g &lt;string&gt; Use only atoms in index group "string" <br>
-n &lt;file&gt; Index file in which to find the group <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: BPR October 2005 


Revision: 1.15-58-g17c5ec7 


## random_box 


Fill a cell with molecules in random orientation using periodic boundary 
conditions. 


Random box can be used to build systems containing a random assortment of 
molecules or to solvate a molecular system. An existing molecule file (solute) 
can also be embedded in the box. The program can also generate molecules 
oriented randomly within a sphere or spherical shell. 


Use: 


### Simple random box 


Generate a periodic cell size 10 x 10 x 10 Ang containing 100 molecules in 
random orientations: 


random_box -x 10 -y 10 -z 10 -n 100 mol.pdb 


### Protein embedded in a box of solvent 


Generate a periodic cell size 100 x 100 x 100 Ang containing a protein 
solvated by water molecules to give the whole system a density of 1kg/dm^3: 


random_box -x 100 -y 100 -z 100 -d 1 -e prot.pdb wat.pdb 


### Mixed random box 


Random_box can handle more than one input molecule. The program will prompt 
for cell size and numbers of each molecule type 


random_box mol1.pdb mol2.pdb 


An alternate way to generate systems containing multiple random components is 
to run random box multiple times, embedding the previous system each time. 


#### Geometric shapes 


Random_box can insert molecules into a variety of geometric shapes; Spheres, 
spherical shells, cylinders and planes. These shapes can be combined by 
specifying multiple shapes in a single run. 


#### Spheres and spherical shells 


Sphere of radius 20 Ang 


random_box mol1.pdb -sphere 20 


Spherical shell of radius 20 Ang and thickness 10 Ang in the centre of a 
periodic cell 


random_box mol1.pdb -sphere 30 -thick 10 


#### Cylinders and rods 


Rod of radius 20 Ang and spanning periodic cell: 


random_box mol1.pdb -cylinder 20 


Cylinder of radius 20 and length 30 Ang 


random_box mol1.pdb -cylinder 20 -thick 30 


#### Planes 


Plane of thickness 30 Ang in the xy plane 


random_box mol1.pdb -plane 30 -xy 


### Holes 


Regularly spaced holes can be inserted into geometric shapes. The number of 
holes and hole radius need to be specified. 


Make a spherical shell with 5 evenly spaced holes with a hole radius of 5 
Angstroms: 


random_box mol1.pdb -sphere 30 -thick 10 -holes 5 -hr 5 


### Multiple Shapes 


Multiple shapes can be inserted into a single periodic cell. The arrangement 
is specified with the -origin flag. The -origin flag can be combined with 
other options to produce multiple more complex structures. 


Currently supported origins are: 


a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A single shape located in the centre of the cell <br>
b&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A single shape located in the centre of the XY plane and at 0 on Z axis <br>
c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Four shapes in a hexagonal pattern in the Z plane <br>

Produce 4 spheres of radius 10 Angstroms with hexagonal packing on the Z plane 


random_box mol1.pdb -sphere 10 -origin c 


### The -invert flag and solvation of shapes 


The -invert flag can be used to select the area not occupied by a shape. This 
is useful for solvation: 


random_box -e sphere_with_holes.pdb -sphere 30 -thick 10 -holes 5 -hr 5 
-invert water.pdb 


### What to do if the required number of molecules won't fit in your desired cell 


If random_box can not fit the desired number of molecules in the cell (i.e. 
the packing process does not complete), try increasing the cell size, followed 
by running an NVT simulation. The system will then shrink to the correct 
volume. 


Flags: 




-e &lt;filename&gt; Embed molecule filename <br>
-n &lt;number&gt; Number of molecules to add <br>
-d &lt;number&gt; Calculate number and add molecules to give this density of entire molecular system (including solute). Note: Number of molecules added is rounded to 3 significant figures. e.g 12345 will be rounded to 12300. <br>
-t&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate embed molecule to centre [origin, cell, none] (default cell) <br>
-rand &lt;number&gt; Randomise input structures by randomly rotating &lt;number&gt; single bonds <br>



Packing flags 


-si &lt;factor&gt; Initial scale factor for atom-atom packing distance (default 0.8) <br>
-sm &lt;factor&gt; Minimum scale factor for atom-atom packing distance (default 0.4) <br>
-incr &lt;number&gt; Number of unsuccesful moleucle packing trials before reducing the scale factor. Default 100 <br>



Geometry Flags 


-x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>
-y&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } Dimensions of the box (Angstroms) <br>
-z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; } <br>
-xy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Build plane or cylinder in the xy plane <br>
-plane &lt;thick&gt; Construct a plane of molecules in the centre of the cell of thicknes &lt;thick&gt;. By default the XZ plane is constructed unless the -xy flag is used. <br>
-sphere &lt;rad&gt; Build sphere with radius &lt;rad&gt; <br>
-cylinder &lt;rad&gt; Build a cylinder with radius &lt;rad&gt; along the X axis <br>
-thick &lt;width&gt; Thickness of spherical shell or length of cylinder (leave unset to generate solid sphere or rod spanning entire periodic cell) <br>
-length &lt;len&gt; Length of cylinder <br>
-holes &lt;number&gt; Number of holes to be evenly distributed on the sphere <br>
-hr &lt;rad&gt; Radius of holes <br>
-origin &lt;type&gt; Origin points for hole vectors: <br>
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; a Single origin in centre of cell (0.5 0.5 0.5) <br>
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; b Single origin at (0.5 0.5 0) <br>
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; c Four origins in xy plane aranged in hexagonal pattern <br>
-invert&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Invert the sense of all 3d structures (i.e. inverting a sphere will give a spherical hole) <br>

Output 


-write&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write intermediate structures <br>
-pc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make files PDB format compatible with maximum SUBID=9999 <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2003 and later 


Revision: 1.76.2.16.2.24 


## Randomise 


Randomise velocities and or coordinates 


Flags: 


-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Randomisation type [vels, coords, vels_coords] (default vels) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.1 


## rdf 


Calculate a radial distribution function 


Flags: 


-i&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (required) RDF Input file <br>

Revision: 1.36.2.1.2.1 


## renumber_residues 


Renumber residues in a file. Number SUBCOUNTs sequentially from 1, and make 
the SUBID for any atom the same as the SUBCOUNT for that atom. Optionally, use 
a different start and a different increment. 


Molecule is sorted before renumbering. All hydrogens are forced to have the 
same residue name, residue number, chain and segid as their parent heavy atom 


Flags: 


-s &lt;number&gt; New starting residue number (default 1) <br>
-i &lt;number&gt; New increment (default 1) <br>
-i &lt;number&gt; New increment (default 1) <br>
-ns&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not sort molecule <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.15-48-g5c3506e 


## residue_distance 


Calculate distance of each residue (or atom) in <comparison_mol> to atoms in 
<source_mol>. Put the data into pdb temperature factor column and or charge 
field 


Useful to calculate how close each protein residue is to a ligand 


Use: residue_distance <source_mol> <comparison_mol1> <comparison_mol2> ... 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Colour by atom rather than residue distance <br>
-k&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Retain only residues (or atoms) closer than 'dist' in output file <br>
-d &lt;dist&gt; Retain dist <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.7.2.1.2.1 


Created: DKC Jan 2008 


## residue_rename 


Change a single residue name 


Rename a single residue type 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename all residues <br>
-e &lt;residue_names&gt; List of residue names to change <br>
-n &lt;numbers&gt; List of residue numbers to change <br>
-r &lt;residue_name&gt; New residue name <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.16.2.1.2.1 


## rms 


Calculate RMS distances between molecules without superimposition 


The first structure in the first file is used as the reference structure. The 
RMS distance is calculated to all subsequent molecules. Output is written to 
<ref_file>.rms 


Heavy atom RMS is calculated by default. The -all flag can be used to include 
hydrogens in the calculations 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use all atoms including hydrogens to calculate RMS. Uses heavy atoms by default <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms into smiles order before doing RMS comparison. This may be useful if molecules have different atom orders. <br>
-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out file containing RMS as SDF_DATA <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2003 


Revision: 1.21.2.1.2.1 


## rms 


Calculate RMS distances between multiple sets of molecules with our without 
superimposition and plot out an RMS matrix 


The -sm flag reorders molecules within the file by name. This is useful if the 
order of molecules is differerent in each file 


The -sa flag sorts atoms within each molecule according to a canonical smiles 
string. This flag deletes all hydrogen atoms and converts all double bonds to 
single bonds. This option allows the calculation of RMS values for different 
tautomers and protonation states. 


Flags: 


-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use all atoms including hydrogens to calculate RMS. Uses heavy atoms by default <br>
-sm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort molecules by name before doing rms comparison <br>
-sa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms into smiles order before doing RMS comparison. This may be useful if molecules have different atom orders. <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2011-14 


Revision: 1.3.2.2.2.1 


## run_mopac 


Minimise a molecule or calculate charges using MOPAC. 


Use: run_mopac [flags] <molecule_file> 


Assumes that the mopac executable is 'mopac', otherwise the executable 
location should be specified on the command line. 


Tested with mopac16 


Default Hamiltonian is: PM7 


Default Mopac options are: NOINTER LET MMOK XYZ GEO-OK. 


Performs a simple check to see if a molecule is a radical. 


To calculate charges on a molecule only use the -1scf flag to set the 1SCF 
option. 


Use the -esp flag to cacluclate ESP charges (not tested lately). 


Flags: 


-ham&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Hamiltonian (default PM7) <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;val&gt; Total charge on molecule <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Radicals are OK. Otherwise, radicals will be skipped <br>
-geo&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use Mopac GEO-OK keyword <br>
-1scf&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; use 1SCF keyword (ie do not minimise) <br>
-ef&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use eigenvector following (EF) <br>
-esp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use esp charges (ESP) <br>
-mozyme&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use MOZYME keyword <br>
-cycles &lt;val&gt; Use CYCLES keyword (default 1000) <br>
-k "string" Additional keywords <br>



-noclean&nbsp;&nbsp;&nbsp;&nbsp; Do not delete Mopac output files <br>
-x&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to Mopac executable <br>



Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 13.6.2000 and on and on 


Revision: 1.1.2.6 


## run_quick 


Minimise a molecule using Quick (https://github.com/merzlab/QUICK) 


Use: runquick [flags] <molecule_file> 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default method is: DFT  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default basis set is: B3LYP  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default calculation is: OPTIMIZE  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Other default options: CUTOFF=1.0d-10 DENSERMS=1.0d-6  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Performs a simple check to see if a molecule is a radical.  <br>

Flags: 


-method&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; HF or DFT (default DFT) <br>
-f &lt;string&gt; Functional [STO-3G, B3LYP, cc-pVDZ, etc]. (Default B3LYP) <br>
-b &lt;string&gt; Basis set [6-311G, 6-311G*, cc-pvDZ, etc] <br>
-c &lt;int&gt; Formal charge <br>
-scf &lt;int&gt; Max cycles for SCF convergence <br>
-densrms &lt;val&gt; DENSRMS option [1.0d-6] <br>
-cutoff &lt;val&gt; CUTOFF option [1.0d-6] <br>
-scf &lt;int&gt; Maximum steps for SCF convergence (Default 30) <br>
-noopt&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not use OPTIMIZE keyword <br>
-steps &lt;int&gt; Steps for OPTIMIZE keyword (0 gives Quick default) <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>

Run options 


-cuda&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use CUDA version <br>
-mpi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use MPI version <br>
-n &lt;val&gt; Number of CPUs to use with MPI version (default 6) <br>
-x &lt;path&gt; Path to Quick executable <br>



Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 17.10.2021 and on and on 


Revision: 1.1.2.8 


## sdf_merge 


Merge data from files based on structure or compound name 


Flags: 


-u&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Merge using this data field (default NAME) <br>
-sum&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sum these fields. Must be suffixed with index (.a, .b, ...). Separated by '|'; <br>
-min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Lowest value from these fields <br>
-s&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort output using this field (default glide_gscore.Rank.best) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 and later 


Revision: 1.1.2.4 


## sdfsplit 


Fast split program to divide an sdf file into smaller files without parsing 
the file 


Use: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sdfsplit file.sdf -s fnumbered produces file_00001.sdf, file_00002.sdf ...  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sdfsplit file.sdf -s fnumbered_2 produces file_01.sdf, file_02.sdf ...  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sdfsplit file.sdf -s numbered produces file_1.sdf, file_2.sdf ...  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sdfsplit -r &lt;CODE&gt; renames molecules using the data in field &lt;CODE&gt;  <br>

Flags: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-n &lt;num&gt; Number of structures to write to each output file (default 1)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-nf &lt;num&gt; Number of files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-f &lt;field&gt; Split structures into files by specified field  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-r &lt;datafield&gt; Rename molecules using the given SDF data field  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-s (numbered / fnumbered) Output style  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-d &lt;dirname&gt; Output directory  <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.20.2.1.2.6 


## Smooth 


Smooth a trajectory by averaging frames without generating 'jumping molecule' 
artefacts caused by molecules crossing the periodic cell boundary 


Flags: 


-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Smoothing window (number of frames, default 2) <br>
-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not average residues that move more than this distance (default 15 Ang) <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC November 2009 


Revision: 1.8.2.2.2.2 


## Starmaker - An all-purpose dendrimer builder 


Starmaker is a command-line script that assembles arbitray dendrimers. 
Monomers are added layer-by-layer starting from a core molecule. The dendrimer 
can be minimised after addition of each monomer layer using the external 
program Macromodel (Schrodinger) 


### Main features: 


The dendrimer is built layer by layer. The first residue in the list is named 
COR. Subsequent residues are named GAA, GAB, GAC, etc. Amide bonds in each 
monomer are recognised and converted to a trans geometry. A repulsion 
potential and random torsional search is used to force monomers to grow in a 
extended geometry. 


At the completion of each layer the molecule is minimised using Macromodel 


### Additional features: 


A file 'stop' in the working directory will terminate the program 


### Use: 


Use: starmaker5 core.mol2 monomer1.mol2 monomer2.mol2 monomer3.mol2 ... 
cap.mol2 


Each mol2 file should contain the monomer for that generation. Each monomer 
needs to have attachment points. These are hydrogen atoms named Q1 (point to 
branch FROM) or Q2 (point to attach TO). The first (core) subunit should 
contain at least one Q1 hydrogen. A standard monomer should contain one Q2 
hydrogen and at least one Q1 hydrogen. A capping group should contain one Q2 
hydrogen. 


Outputs a mol2 file 'final_dendrimer.mol2' and intermediate states as 
layer_XX.mol2 


### Flags 


-rc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force first residue to be renamed to 'COR'. <br>
-clash&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distance below which two atoms clash <br>
-noh&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not include hydrogens in geometry optimisation (default) <br>
-cut&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Distance above which the repulsive potential is ignored <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; External minimisation prog (Macromodel) <br>
-iter&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of iterations when optimising geometry (no clashes) <br>
-conv&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of steps with the same best coordinates before optimisation stops <br>
-min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of steps <br>
-wc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out structure after each monomer to file current.mol2 <br>
-wi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write out structure after every optimisation step <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.15-89-gc2af617 


## Tabulate 


Output data from molecule file (e.g. .sdf or .maegz file) and create delimited 
text file. Optionally plot histograms of data values using gnuplot 




Flags: 


Created: DKC 2003 and on 


Revision: 1.31.2.2.2.7 


## vesicle_builder 


Build a vesicle using periodic boundary conditions. 


Lipid molecules are reoriented along their polar axis, unless the -no flag is 
specified 


Usage: vesicle_builder <lipid_file> -x <cellx> -y <celly> -z <cellz> -n1 
<number of lipids layer1> [ -n2 <number of lipids layer2> -e <molecule to 
embed> ] 


Flags: 


-n1 &lt;number&gt; Number of instances of lipid to use to create monolayer 1 (outer) <br>
-n2 &lt;number&gt; Number of instances of lipid to use to create monolayer 2 (inner) <br>
-e &lt;filename&gt; Embed molecule filename <br>

Geometry flags 


-rad &lt;dist&gt; Outside radius of vesicle <br>
-x &lt;dist&gt; | <br>
-y &lt;dist&gt; | Dimensions of the cell in each direction <br>
-z &lt;dist&gt; | <br>
-no&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not reorient the polar axis of embedded lipid molecules <br>

Packing flags 


-rand &lt;number&gt; Number of lipid bonds to randomise <br>
-trans&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Translate the embedded molecule to the cell centre <br>
-es &lt;number&gt; Additional distance addedd between the mid-points of the bilayers <br>
-si &lt;factor&gt; Initial scale factor for atom-atom packing distance (default 0.8) <br>
-sm &lt;factor&gt; Minimum scale factor for atom-atom packing distance (default 0.4) <br>
-incr &lt;number&gt; Number of unsuccesful moleucle packing trials before reducing the scale factor. Default 100 <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 and later 


Revision: 1.1.2.9 


## water_to_ion 


Replace random water molecules (single atom) counter ions. Defaults to Na 
ions. Covers a very large range of inorganic ions. 


Note: sorting atoms (the default) can cause problems 


Flags: 


-e &lt;string&gt; Ion element (default Na). Now handles a large range of single-atom cations and anions. <br>
-n &lt;number&gt; Number of ions to add <br>
-r &lt;string&gt; Residue name of added ions (optional) <br>
-a &lt;string&gt; Atom name to use for ions (optional) <br>
-amber&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add PDB TER record after each ion to produce input for Amber <br>
-sort&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort the atom order in the resulting output file <br>
-ra&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Renumber atoms <br>
-rr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Renumber residues <br>
-rand&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Replace random water molecules distributed through water molecules in file. Otherwise the first 'n' waters in the file are selected. <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: David K. Chalmers & Benjamin P. Roberts 2000-8 


Revision: 1.42.2.1.2.2 


## write_g03zmat 


Read a molecule and turn it into a Z-matrix, with atomic connectivity given in 
the order in which the atoms appear. 


Usage: write_g03zmat <file1> <file2> ... 


Created: Benjamin P. Roberts, October 2008 


Revision: 1.1.2.2.2.1 


## write_gro 


Read any silico format and write a xyz file. 


Use: write_gro <file1> <file2> ... 


Flags: 


-split&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split multiple molecule file into individual structures <br>
-prec&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Precise - increase precision of .gro file to 0.001 Angstroms <br>
-vprec&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Very precise - increase precision of .gro file to 0.0001 Angstroms <br>
-zv&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Zero velocities - set all velocities to zero <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2.2.5 


Created: DKC November 1999 and later. 


## write_itp 


Read and Write Gromacs .itp topology files 


Flags: 


-coord&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;file&gt; Add coordinates from &lt;file&gt; into itp molecule. The coordinate file and ITP file are checked for consistency. <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Replace peptide backbone parameters with standard gromacs types <br>
-heavy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Multiply hydrogen mass by a factor of 4 and subtract the additional mass from the parent heavy atom (make heavy hydrogens) <br>
-resname&nbsp;&nbsp;&nbsp;&nbsp; Set residue name <br>
-sort&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms into canonical SMILES order (requires coord file) <br>
-rename&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename atoms using a simple naming scheme (requires coord file) <br>
-rmethod&nbsp;&nbsp;&nbsp;&nbsp; Renaming method (see mol_rename_atoms() subroutine) [consecutive, consecutinve-h, element, element-h', 1, 'element-h <br>
-dump&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dump all bonds, angles, dihedrals, etc as individual molecules to itp-data directory to a assist with debugging toplogies <br>

Created: DKC 2001 and on... 


Revision: 1.1.2.2.2.16 


## write_mmod 


Read a molecule and write a Macromodel format file. 


Useage: write_mmod <file1> <file2> .... 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 and on 


Revision: 1.1.2.1 


## write_merck 


Read a molecule file and write a Merck format file 


Use: write_merck <file1> <file2> ... 


Flags: 


-b,--regenerate-bondorders Regenerate bondorders <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 


Revision: 1.1.2.2.2.1 


## write_mmod 


Read a molecule and write a Macromodel format file. 


Useage: write_mmod <file1> <file2> .... 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 


Revision: 1.1.2.2.2.1 


## write_mol 


Read a molecule and write (by default) in the same format. 


Flags: 


-check&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Check integrity of molecules <br>
-merge&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Merge all molecules into a single molecule <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2001 


Revision: 1.30.2.1.2.4 


## write_mol2 


Read any silico format and write a mol2 file 


## Use 


write_mol2 <file1> <file2> ... 


Flags: 


-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename molecules using SDF data field for conversion from SDF to mol2 <br>
-mm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write output in MolMol Mol2 format. <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write output in Mol2 Protein format. <br>
-dr&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Print additional debugging information for rings <br>
--nosingle Do not use single molecule read/write routines <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.15-21-g2cee521 


Created: DKC November 1999 and later. 


## write_mopac 


Read any Silico format and write a MOPAC Cartesian file. 


Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000. 


Revision: 1.1.2.2.2.1 


## write_pdb 


Read any Silico format and write a pdb file 


Flags: 


-de &lt;alternate_conf&gt; Delete all disorderd (ALT) atoms except those in the specified alternate conformation (default Conformation 'A') <br>
-dl &lt;value&gt; Delete low occupancy atoms. Atoms with occupancy below &lt;value&gt; are deleted <br>
-model&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write each structure as a separate MODEL <br>
--nosingle Do not use single molecule read/write routines <br>
-debug&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC November 1999-2025 


Revision: 1.1.2.2.2.9 


## write_sdf 


Read any silico format and write an sdf file. 


Optionally add a 'name' field, rename the molecule or remove SDF data 


Flags: 


-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Number of structures to read <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; SDF data field to use if renaming molecules <br>
-a&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add 'NAME' field to SDF_DATA using molecule name encoded in the first line of the file <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Add 'FILENAME' field to SDF_DATA using using input filename <br>
-clean&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Remove all sdf data (except name) <br>
-noparse&nbsp;&nbsp;&nbsp;&nbsp; Do not parse SDF data (Only works with SDF input files) <br>
-nosingle&nbsp;&nbsp; Do not use 'single' molecule reading routines <br>
-2D&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force structure to be written as 2D <br>
-3D&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force structure to be written as 3D <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 - 2002 


Revision: 1.1.2.2.2.6 


## write_seq 


Sequence format test script 


Silico protein/DNA sequence format routines are under development and 
incomplete 


Flags: 


-c,--combine Combine sequences from all files to a single file <br>
-n,--number-of-residues Number of residues per line in output <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.30.2.2.2.2 


## write_tab 


Read and write character delimited files (default space delimited) 


Flags: 


-d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Delimiter (tab [\t], space [\s+], comma [,]) <br>
-smooth&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Smooth data <br>
-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Smoothing window <br>
-noheader&nbsp;&nbsp; First line is NOT header line <br>

Created: DKC November 2003 


Revision: 1.1.2.2.2.2 


## write_tinker 


Read any Silico format and write a tinker xyz file. 


Usage: write_tinker <file1> <file2> ... 


Flags: 


-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force field file (currently supports opls) <br>
-n&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use normal force field atom types (not special Tinker types) <br>
-p&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Location of master parameter database <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename atoms so atom name = atom type <br>
-w&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write info on unassigned atoms to standard output <br>
-v&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Verbose output, equivalent to -r -w <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Created: DKC 2000 


Revision: 1.1.2.2.2.2 


## write_top 


Read and Write Gromacs .top and .itp topology files 


Flags: 


-coord&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &lt;file&gt; Add coordinates from &lt;file&gt; into itp molecule. The coordinate file and ITP file are checked for consistency. <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Replace peptide backbone parameters with standard gromacs types <br>
-heavy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Multiply hydrogen mass by a factor of 4 and subtract the additional mass from the parent heavy atom (make heavy hydrogens) <br>
-resname&nbsp;&nbsp;&nbsp;&nbsp; Set residue name <br>
-sort&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Sort atoms into canonical SMILES order (requires coord file) <br>
-rename&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Rename atoms using a simple naming scheme (requires coord file) <br>
-rmethod&nbsp;&nbsp;&nbsp;&nbsp; Renaming method (see mol_rename_atoms() subroutine) [consecutive, consecutinve-h, element, element-h', 1, 'element-h <br>
-dump&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dump all bonds, angles, dihedrals, etc as individual molecules to itp-data directory to a assist with debugging toplogies <br>

Created: DKC 2001 and on... 


Revision: 1.1.2.1 


## write_xyz 


Read any silico format and write an xyz file. 


Use: write_xyz <file1> <file2> ... 


Flags: 


-split&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Split multiple molecule file into individual structures <br>

Some general Silico Flags are -o, -O, -compress, -nosort, -notype, -fast, 
-bond, -nobond, -noconect, -ms, -me, -ss, --last-structure, -longres --timing, 
--quiet, --very-quiet and -debug. These are available in most scripts. Some 
special options apply to SDF and PDB files. 


For more information rerun this program using '--help-flags' 


Revision: 1.1.2.2 


Created: DKC November 1999 and later. 


## Directory: /vcp1/people/david/bin/Silico1.15/data 


## amino_acid_atoms.dat 


Amino acid connectivities and atom names 


## amino_acids_smiles.dat 


Data file containing SMILES strings corresponding to amino acid residues. The 
following conventions are used for protonation-state suffixes: 


Side chains ("S" column) 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "0": standard.  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "S": Substitution on terminal S, e.g., in disulfide (CYS).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "X": Deprotonation on NZ (LYS).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "F": Deprotonation on NE (ARG).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "J": Deprotonation on NH (ARG).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "D": Protonation on ND but not NE (HIS).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "E": Protonation on NE but not ND (HIS).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "H": Protonation of the side-chain carboxylate (ASP, GLU).  <br>

Backbones ("B" column) 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "0": standard (intrapeptide).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "N": N-terminal amino acid, protonated on terminal n  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "T": n-terminal amino acid, neuTral on terminal n  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "C": C-terminal amino acid, deprotonated on terminal carboxylate  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "O": c-terminal amino acid, neutral  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "A": Acidic isolated amino acid (fully protonated)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "B": Basic isolated amino acid (fully deprotonated)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "U": isolated amino acid with neUtral backbone  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- "Z": isolated amino acid with Zwitterionic backbone  <br>

Atom names are given next. These are the standard atom names, as obtained from 
the PDB specifications. The following variations are observed (where the PDB 
is not specific enough): 


- Arginine: where deprotonation has occurred on an eta nitrogen, the 
deprotonated nitrogen is called "NH1" (since it has 1 proton attached). Its 
counterpart with two protons is called "NH2". 


Data file containing SMILES strings for amino acid residues. 




#### Atom specifiers 


-as&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Atom specifier strings (multiple -as flags are allowed) <br>

Several flags are available as atom specifier shortcuts 


-back&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Backbone atoms <br>
-ca&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA atoms <br>
-cacb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CA and CB atoms <br>
-heavy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Nonhydrogen atoms <br>
-sas&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Short specifier [ca, cacb, back, back_disulfide, heavy] <br>
-f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Read atom specifiers from a file <br>

#### Atom specifier examples: 


ANAME:CA&nbsp;&nbsp;&nbsp;&nbsp; All atoms called 'CA'. <br>
ANAME:CA,CB,CG Returns all atoms called 'CA' or 'CB' or 'CD'. <br>
ANAME:CA,RESNAME:TRP All atoms called 'CA' in all residues called TRP <br>
ANAME:CA,SUBID:4 All atoms called 'CA' in residue number 4 <br>
ELEMENT:!H All nonhydrogen atoms <br>
SEGID:PROT,SEGID:LIG All atoms with the SEGID set to PROT or LIG <br>

Successive atom specifications can be made. Each is separated by a '|' 


ANAME:CB|ANAME:CA|ANAME:CD Returns all atoms called 'CB', 'CA', 'CD'. <br>
ANAME:CA,RESNUM:1|ANAME:CA,RESNUM:4 Returns 'CA' atoms from residue 1 and 4. <br>

Atom specifiers are case sensitive. 


## charmm_amino_acid_atoms.dat 


Amino acid connectivities and names according to Charmm22 


Revision: 1.7.2.1.2.1 


## cns_amino_acid_atoms.dat 


Connection tables and atom names of amino acids according to 'CNS' 


## cyana_amino_acid_atoms.dat 


Connection tables and atom names of amino acids according to Cyana2 


Revision: 1.7.2.1.2.1 


Standard fragment flags 


-fp _max&nbsp;&nbsp;&nbsp;&nbsp; Max fragment size <br>
-fp _min&nbsp;&nbsp;&nbsp;&nbsp; Min fragment size <br>
-fp _type&nbsp;&nbsp; Fragment typing method (simple (default), steric) <br>
-fp _h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include hydrogens in fragments (none (default), polar, all) <br>

The list of standard gromacs flags 


Gromacs settings: 


-g _exe_&lt;string&gt; Gromacs executable (default gmx) <br>
-g _nt_&lt;int&gt; Gromacs -nt flag <br>
-g _maxwarn_&lt;num&gt; Gromacs grompp maximum number of warnings <br>
-g _pin&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Gromacs pin threads flag <br>
-g _pinstride_&lt;int&gt; Gromacs pinstride for use with gmx-pin <br>
-g _pinoffset_&lt;arr&gt; Space delimited array of pinoffset values <br>
-g _f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Additional gromacs flags <br>

Gromacs local job options: 


-h &lt;hostlist&gt; Host list for local jobs. Space separated string. Use quotes on the command line. <br>
-maxj &lt;int&gt; Maximum number of GPU jobs(gmx, desmond, etc) allowed per GPU <br>

SLURM job options: 


-s _qos_&lt;string&gt; SLURM quality of services <br>
-s _time_&lt;string&gt; SLURM maximum job run time (hrs:mins:secs) <br>
-s _nnode_&lt;num&gt; SLURM number of nodes to request (default 1) <br>
-s _ncpu_&lt;num&gt; SLURM number cpus <br>
-s _ngpu_&lt;num&gt; SLURM number of GPUs <br>
-s _part_&lt;string&gt; SLURM partition <br>

The list of Mopac flags 


-ham&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Method [AM1, PM3, etc] <br>
-c&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Formal charge <br>
-r&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Run calculations on radical species <br>
-geo&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Run calculations on radical species <br>
-1scfCalculate charge only (1SCF) <br>
-ef&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Geometry optimisation using eigenvector following <br>
-esp&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use ESP charges <br>
-mozyme&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Use mozyme keyword <br>
-gnorm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Value for GNORM keyword <br>
-k&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Additional keywords <br>
-exe&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Location of MOPAC executable <br>

## opls_amino_acid_atoms.dat 


Data file containing mappings of amino acid atoms to OPLS-AA atom types. This 
file assumes that atoms have been named according to the standards in the file 
amino_acids_smiles.dat (derived in turn from the PDB standard). 


Not all amino acid atoms have been included here; only those with OPLS-AA 
types that would not otherwise be assigned to them. 


Hydrogens are not included. 


## opls_hydrogens.dat 


Data file containing mappings of OPLS-AA hydrogen atom types to their parent 
heavy atoms. 


## oplsaa_parameters3.dat 


Silico format parameter file for the OPLS-AA force field (version 3) 


This file is used by the OPLS atom typing and psf/prm generation routines. 


Note: The file is tab delimited. Any altered lines must use tabs to separate 
fields. 


The ordering of dihedrals in the file is important. Dihedrals are typed using 
the first matching set. Therefore specific types should be listed first and 
more general matches should be placed later in the file. Dihedrals with an "@" 
at the start of the line are duplicate entries, which are to be reported if 
dihedrals in the molecule match whatever they are duplicates of. This happens 
so potential misassignments are drawn to the user's attention and he has the 
opportunity to do something about them. 


Most of these parameters were obtained from W. L. Jorgensen in 2004. These are 
noted in the source column as "Jorgensen OPLS-AA Parameters, June 2004". 




Version 3 


This now includes specific atom types in dihedral types. Eg dihedrals can now 
be defined as CT-CT-2185-CT. 


Carbohydrate atom types have been added using additional atom types to 
differentiate carbohydrates from other polyols. 


OPLS-AA/L dihedral types have been added for polypeptides.  4 Primary 
References 


1. W. L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, "Development and 
Testing of the OPLS All-Atom Force Field on Conformational Energetics and 
Properties of Organic Liquids", J. Am. Chem. Soc., 117, 11225-11236 (1996) 


2. D. S. Maxwell, J. Tirado-Rives and W. L. Jorgensen, "A Comprehensive Study 
of the Rotational Energy Profiles of Organic Systems by Ab Initio MO Theory, 
Forming a Basis for Peptide Torsional Parameters", J. Comput. Chem, 16, 
984-1010 (1995) 


3. W. L. Jorgensen and N. A. McDonald, "Development of an All-Atom Force Field 
for Heterocycles. Properties of Liquid Pyridine and Diazenes", THEOCHEM-J. 
Mol. Struct., 424, 145-155 (1998) 


4. N. A. McDonald and W. L. Jorgensen, "Development of an All-Atom Force Field 
for Heterocycles. Properties of Liquid Pyrrole, Furan, Diazoles, and 
Oxazoles", J. Phys. Chem. B, 102, 8049-8059 (1998) 


5. R. C. Rizzo and W. L. Jorgensen, "OPLS All-Atom Model for Amines: 
Resolution of the Amine Hydration Problem", J. Am. Chem. Soc., 121, 4827-4836 
(1999) 


6. M. L. P. Price, D. Ostrovsky and W. L. Jorgensen, "Gas-Phase and 
Liquid-State Properties of Esters, Nitriles, and Nitro Compounds with the 
OPLS-AA Force Field", J. Comput. Chem., 22, 1340-1352 (2001)  4 Secondary 
References 


7. K. Khan and T. C. Bruice, "alpha-Ketoamides and alpha-Ketocarbonyls: 
Conformational Analysis and Development of All-Atom OPLS Force Field", Bioorg. 
Med. Chem., 8, 1881-1891 (2000) 


8. W. Damm, A. Frontera, J. Tirado-Rives and W. L. Jorgensen, "OPLS All-atom 
force field for carbohydrates. J. Comp. Chem. 18, 1955-1970 (1997) 


9. W. L. Jorgensen and D. L. Severance, "Aromatic-Aromatic Interactions: Free 
Energy Profiles for the Benzene Dimer in Water, Chloroform and Liquid 
Benzene", J. Am. Chem. Soc. 112, 4768-4774 (1990) 


10. W. L. Jorgensen, J. D. Madura and C. J. Swenson, "Optimized Intermolecular 
Potential Functions for Liquid Hydrocarbons", J. Am. Chem. Soc. 106, 6638-6646 
(1984) 


11. J. M. Briggs, T. B. Nguyen and W. L. Jorgensen, "Monte Carlo Simulations 
of Liquid Acetic Acid and Methyl Acetate with the OPLS Potential Functions", 
J. Phys. Chem. 95, 3315-3322 (1991) 


12. W. L. Jorgensen, "Optimized Intermolecular Potential Functions for Liquid 
Alcohols", J. Phys. Chem. 90, 1276-1284 (1986) 


13. W. L. Jorgensen, "Intermolecular Potential Functions and Monte Carlo 
Simulations for Liquid Sulfur Compounds", J. Phys. Chem. 90, 6379-6388 (1986) 


#### Sphere packings  Source; 


1 .. 13 Optimal packings of a sphere > 130 Icosahedral packings  N. J. A. 
Sloane Email address: njasloane@gmail.com Home page: NeilSloane.com/ 


## pdb_space_groups.dat 


Data file containing mappings of Silico space groups to their corresponding 
PDB format strings. 


The list of standard Silico flags (i.e., those found in most Silico programs). 


Standard Silico Flags: 


-o &lt;format&gt; Output file format <br>
-O &lt;file&gt; Output file name <br>
-nosort&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not sort atoms before output <br>
-notype&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not type atoms before output <br>
-fast&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fast write (equivalent to -nosort -notype) <br>
-bond&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force regeneration of bond records <br>
-nobond&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not regenerate bond records <br>
-noconect&nbsp;&nbsp; Do not generate CONECT records for PDB output <br>
-compress &lt;method&gt; Output file compression. Method = gz or bz2 <br>
--timing&nbsp;&nbsp;&nbsp;&nbsp; Print more detailed timing information <br>
--quiet&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Quiet (suppress notes) <br>
--very-quiet Very quiet (suppress notes and warnings) <br>
-debug&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write debugging messages   <br>

General flags implemented in some programs: 


-me&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum energy relative to first structure <br>
-ms&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of structures to read <br>
-ss&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Start reading from this structure <br>
--stride &lt;n&gt; Read every nth structure from file <br>
--unique-atom-names Make atom names in molecules unique  <br>

PDB file specific flags: 


-dupconect Represent multiple bonds using duplicated CONECT records in PDB output <br>
-longres&nbsp;&nbsp;&nbsp;&nbsp; Use 4 character residue names in pdb files  <br>

SD file specific flags: 


--sdf-2D&nbsp;&nbsp;&nbsp;&nbsp; Force SDF files to 2D format <br>
--sdf-3D&nbsp;&nbsp;&nbsp;&nbsp; Force SDF files to 3D format <br>
--sdf-squash-hcount Do not write any values to the SDF hcount field  <br>

Gaussian and Gamess specific flags: 


--last-structure Retain only the last read structure (implemented in read_gaussian and read_gamess) <br>

The list of standard Silico flags (i.e., those found in most Silico programs). 


Standard Silico Flags: 


-o &lt;format&gt; Output file format <br>
-O &lt;file&gt; Output file name <br>
-nosort&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not sort atoms before output <br>
-notype&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not type atoms before output <br>
-fast&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Fast write (equivalent to -nosort -notype) <br>
-bond&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Force regeneration of bond records <br>
-nobond&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Do not regenerate bond records <br>
-noconect&nbsp;&nbsp; Do not generate CONECT records for PDB output <br>
-compress &lt;method&gt; Output file compression. Method = gz or bz2 <br>
--timing&nbsp;&nbsp;&nbsp;&nbsp; Print more detailed timing information <br>
--quiet&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Quiet (suppress notes) <br>
--very-quiet Very quiet (suppress notes and warnings) <br>
-debug&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Write debugging messages   <br>

General flags implemented in some programs: 


-me&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum energy relative to first structure <br>
-ms&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of structures to read <br>
-ss&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Start reading from this structure <br>
--stride &lt;n&gt; Read every nth structure from file <br>
--unique-atom-names Make atom names in molecules unique  <br>

PDB file specific flags: 


-dupconect Represent multiple bonds using duplicated CONECT records in PDB output <br>
-longres&nbsp;&nbsp;&nbsp;&nbsp; Use 4 character residue names in pdb files  <br>

SD file specific flags: 


--sdf-2D&nbsp;&nbsp;&nbsp;&nbsp; Force SDF files to 2D format <br>
--sdf-3D&nbsp;&nbsp;&nbsp;&nbsp; Force SDF files to 3D format <br>
--sdf-squash-hcount Do not write any values to the SDF hcount field  <br>

Gaussian and Gamess specific flags: 


--last-structure Retain only the last read structure (implemented in read_gaussian and read_gamess) <br>

## Directory: /vcp1/people/david/bin/Silico1.15/lib 


## silico.pm 


The parent silico module 




Silico: a perl molecular toolkit. 




Written by David Chalmers. 29.1.2000 and onward. 




Silico.pm is the Silico parent module. 


Loads a library of of subroutines for handling molecular structures and 
routines to read and write a number of standard molecular modelling file 
formats. 




Silico.pm is usually loaded by the subroutine silico_setup which defines the 
location of the Silico libraries and executables. 




By default the following Silico modules are loaded: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_defnitions - info about atoms etc  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_check - molecule checking routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_control - program control and user interface routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_data - internal data handling routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_debug - debugging routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_doc - automatic documentaton generation routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_geom - molecular geometry routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_hydrogens - Add and delete hydrogens from molecules  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_io - general input/output routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_mol2 - read/write Sybyl mol2 files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_molecules - general operations on molecules  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_pdb - read/write Brookhaven pdb files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_prop - calculate molecular and atomic properties  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_rings - identify rings in molecules  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_sdf - read/write MDL sdf files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_statistics - statistical routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_residue - residue routines  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_asl - atom selection language   <br>

These are loaded by silico_io if required 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_gromacs - Gromacs format files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_merck - Merck format files  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_tinker - Tinker format files  <br>



These have to be specified using a 'require' statement 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_smiles - Read pseudo smiles and type atoms  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_sequence - Extract protein sequences  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;silico_split - Split up molecules based on connectivities   <br>

The global variable $silico_version is also set. 


Revision: 1.27.2.2.2.6 


## silico_asl.pm 


Silico atom selection routines (atom specifier) 


Revision: 1.15-39-g4f8dfe8 


## silico_check.pm 


Silico routines to check the quality of molecules 


Revision: 1.34.2.1.2.3 


## silico_control.pm 


Various subroutines for program control, user interface, etc. 


Revision: 1.15-75-gd957e88 


## silico_data.pm 


Various subroutines for handling text data and other data structures. 


Revision: 1.18.2.1.2.5 


## silico_debug.pm 


Silico general debugging routines. 


Revision: 1.35.2.3.2.4 


## silico_definitions.pm 


Define Silico global constants. 


This module contains no subroutines but defines severable global arrays. These 
are: 




@Atomic_Elements - an array containing the periodic table. 


@Atomic_Masses - array containing atomic masses. 


@Amino_Acids - the standard amino acids. 


@Amino_Acids - the standard amino acids. 


@Amino_Acids1 - the standard amino acids 


@Nucleic_Acids - the standard nucleic acids. 


@Nucleic_Acids1 - the standard nucleic acids - lower case one letter code. 


%Amino_Nucleic_Acids - hash with AA and NA three-letter code as key with 
one-letter codes as value. 


%Amino_Nucleic_Acids1 - hash with AA and NA one-letter code as key with 
three-letter codes as value. 


@Vdw_Radii - VdW radii for determinining if two atoms are bonded 


Revision: 1.32.2.3.2.1 


## silico_doc.pm 


Silico routines to print formatted documentation in text and html formats 
which is generated from program comments that have been marked up using a 
simple syntax 




Simplified markup syntax: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#&lt; = Start marked up comments  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#&gt; = End marked up comments  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#! = Heading/Title (Heading 2)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#!1, #!2, #!3, #!4 Different heading levels  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#? Short description of program/subroutine used for indexing.  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#. = New paragraph.  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#, = Indent  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#; = Hanging indent (long).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#: = Hanging indent (short, for numbering).  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#- = List.  <br>



Standard headings: 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#F = Flags  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#SF = Silico Flags  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;#GF = Gromacs Flags  <br>



Revision: 1.36.2.4.2.7 


## silico_fp.pm 


Routines to generate molecular fragments and make Tanimoto comparisons etc 


Revision: 1.46.2.4.2.1 


## silico_gamess.pm 


Silico routines for reading and writing Gamess files 


Revision: 1.7.2.1.2.1 


## silico_gaussian.pm 


Silico routines for reading and writing Gaussian Z-matrix and Cartesian files 


Revision: 1.8.2.1.2.2 


## silico_geom.pm 


Silico routines to calculate or operate on molecular geometries (distances, 
angles, etc) 


Revision: 1.15-66-g8b88391 


## silico_gnuplot.pm 


Silico wrapper routines to run gnuplot 


Revision: 1.1.2.2.2.1 


## Silico_gromacs.pm 


Routines specific to Gromacs 


Revision: 1.83.2.10.2.50 


## silico_hb.pm 


Silico hydrogen bond routines 


Revision: 1.1.2.2 


## silico_hydrogens.pm 


Routines to add and delete hydrogen atoms 


Revision: 1.46.2.7.2.9 


## silico_io.pm 


Silico general file handling routines and routines for getting arguments and 
flags from the command line. 


Revision: 1.15-74-g1e475b1 


## silico_merck.pm 


Routines to read and write Merck format files for use with MMFF94 in Charmm 


Revision: 1.34.2.1.2.2 


## silico_mmod.pm 


Silico routines specific to macromodel format files 


Revision: 1.15-73-g3bc62ee 


## silico_mol2.pm 


Silico routines to deal with Sybyl mol2 files. 


Revision: 1.15-46-gd5da50c 


## silico_molecules.pm 


Silico routines which operate on any atom, molecule, or ensemble. 


Revision: 1.15-67-gf7b58c9 


## silico_mopac.pm 


Silico Mopac and Divcon routines. 


Revision: 1.56.2.1.2.11 


## silico_pdb.pm 


Silico routines to read and write pdb format files. 


Revision: 1.15-57-gb152f62 


## silico_prop.pm 


Silico routines to calculate molecular and atomic properties including 
molecular interface interactions 


Revision: 1.15-72-gda9b5d8 


## silico_quick.pm 


Silico Quick routines. 


Revision: 1.1.2.11 


## silico_rama.pm 


Silico routines to calculate Ramachandran plots and peptide dihedral angles 


Revision: 1.12.2.4.2.7 


## silico_rbox.pm 


Subroutines for random_box and bilayer_builder 


Revision: 1.1.2.17 


## silico_residue.pm 


Silico routines for handling residues in molecules 


Revision: 1.1.2.11 


## silico_rgmx.pm 


Silico routines to run Gromacs 


Revision: 1.1.2.33 


## silico_rgmx.pm 


Silico routines to run Gromacs 


Revision: 1.1.2.13 


## silico_rings.pm 


Silico routines to identify rings in molecules. 


Revision: 1.15-42-g45ba5cf 


## silico_schedule.pm 


Silico routines to schedule jobs on local workstations 


Revision: 1.1.2.31 


## silico_sdf.pm 


Silico routines to read and write MDL sdf and rdf files. 


Revision: 1.69.2.8.2.8 


## silico_sdf.pm 


Silico routines to read and write MDL sdf and rdf files. 


Revision: 1.69.2.8.2.8 


## silico_sequence.pm 


Routines to deal with protein/dna sequences 


Sequences are stored as hashes in the same manner as Silico molecules. 


Sequence details 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{NAME} sequence name  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{NUMRESGAP} number of residues in sequence including gaps  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{NUMRES} number of residues in sequence excluding gaps  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{SEQ}[] array of residues  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{ENS} associated 3D molecules (pointer to ensemble)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{TYPE} type of information stored in sequence record  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$seq-&gt;{SOURCE} format of data source   <br>

Standard residue details ($seq->{TYPE} = 'STD') 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA1} One-letter code  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA3} Three-letter code  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{NUM} Residue number  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{ATOMS} 3D residue (array of atoms)   <br>

Sequence with associated structure ($seq->{TYPE} = 'STRUCT') contains standard 
residue details plus 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{ATOMS} Atoms associated with residue   <br>

Consensus residue details (seq->{TYPE} = 'CONSENSUS') 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA1} One-letter code (Consensus sequence)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA3} Three-letter code (Consensus sequence)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{NUM} Residue number  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{ATOMS} 3D residue (array of atoms)  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AAS}[] Array of residues found at this position  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA3S}[] Array of residues found at this position  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{FREQS}[] Frequency of each residue   <br>

Clustal conservation (seq->{TYPE} = 'CLUSTAL_CONS') 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA1} Conservation symbyl (*, : , ., ' ' )  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{AA3} Three-letter code  <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$res-&gt;{NUM} Residue number  <br>

Revision: 1.44.2.1.2.2 


## silico_smiles.pm 


Silico library dealing with SMILES strings 


Revision: 1.20.2.9.2.20 


## silico_solvate.pm 


Routines for solvation, density etc 


Revision: 1.6.2.2.2.7 


## silico_split.pm 


Routines to split a molecule into an ensemble of separate molecules. 


Revision: 1.52.2.5.2.24 


## silico_statistics.pm 


Silico routines to calculate statistical quantities 


Revision: 1.9.2.3.2.1 


## silico_tabledata.pm 


Handle data that is in tabular form 


Revision: 1.33.2.4.2.6 


## silico_tinker.pm 


Silico routines specific to Jay Ponder's TINKER and SLEUTH software 


Revision: 1.40.2.1.2.4 


## silico_xyz.pm 


Silico routines specific to Jay Ponder's TINKER and SLEUTH software 


Revision: 1.1.2.3 

