This is a README file for GECKOV

USAGE:
Launch the GECKOV by typing cd'ing into the directory and entering on the command line:

python GECKOV.py

GECKOV is a software tool designed to create phylogenetic trees using a variety of methods which do not require a multiple sequence alignment (MSA).  Trees are visualized using the open source Archaeopteryx software tool written by Christian M. Zmasek from phylosoft.org.


Add Sequences:
	This button launches a dialog for the user to first select a FASTA file and then computes distances between the sequences for all of the calculation types.  If the number of sequences is large, this may take a while (a few hours).  The calculations are stored in a database, so that they can be quickly retrieved for future use.  This is useful if the user wants to experiment with different calculations or drawing different tree types on the same dataset.

Accession Numbers:
	Accession numbers should be unique.  For each entry in the FASTA file, the accession number or string containing the accession number should immediately follow the '>'. Such as:

>gi|12175745|ref|NC_002645.1| Human coronavirus 229E, complete genome

The sequence name should be separated from the accession number by a space.  

If the accession string contains pipes ('|') like the entry above, then the string will be parsed and the number between the third and fourth pipe will be taken as the accession.  This is the format of the NCBI viral genome database entries and is the preferred format for GECKOV.  If there are not enough pipe characters or the accession number is not between the appropriate pipes, GECKOV may use the wrong number for the accession or return an error.

If the virus in question does not have an NCBI accession number, the user should provide one of their own in a string without pipes such as:

>123456.7 Unknown Coronavirus

If any pipes are in the string, GECKOV will attempt to parse it as an NCBI accession number and an error could result.


Quick Calculation:
	If the user already knows what type of calculation they want to use, the Quick Calculation option is available.  This will take a FASTA file as its input and draw a phylogenetic tree based upon the selected options in the menu on the right.  The OUTGROUP for the tree should be the first sequence in the file.  The calculations are not stored, so every time the user wants to change a parameter, the distances must be recalculated.  Nevertheless, this method is great for the impatient user or for use on large sets of sequences.


Sequence type:
	Currently, all sequences inputted to the GECKOV must be nucleotide sequences.  There are two other sequence types on which GECKOV can operate.  Nucleotide sequences can be translated to amino acid sequences.  Nucleotide sequences can also be generalized into purine bases and pyrimidine bases.  So all, adenines and guanines are translated to "R's" for purines and thymines and cytosines to "Y's" for pyrimidines.  Select the appropriate radio button of the type of sequence you want use.

Calculation Type:
	Frequency Based:
		Odds Ratio:  The odds ratio is a measure of bias between the observed and expected k-nucleotide frequencies: r_odds = (freq_obs/freq_exp)
		Log Odds: The logarithm of the odds ratio: r_logodds = log(freq_obs/freq_exp)
		Poisson Deviates: Measure of bias from the normalized deviation from a Poisson distribution: r_poisson = (freq_obs - freq_exp)/(sqrt(n*freq_exp*(1-freq_exp))), where n is the total observed k-nucleotides
		Odds Difference: Similar to Poisson, but with the normalization factor: r_oddsdiff = freq_obs - freq_exp
		Jensen-Shannon: Creates a frequency profile of k-nucleotides in the sequence.  Uses the Kullback-Leibler divergence to calculate the Jensen-Shannon divergence as a measure of the distance between two probability distributions
	Vector Based:
		Moment Vector:
			Moments can be used to characterize a graphical curve.  A graphical representation of the DNA or amino acid sequence based on a vector system in which each of the four nucleotides corresponds to a different vector in the first or fourth quadrant of the Cartesian coordinate system.  Each sequence has a N-dimensional moment vector associated with it, but two or three moments will characterize the protein sufficiently well.  The distance between two moment vectors can be calculated and a phylogenetic tree created. 
		Natural Vector:  The natural vector is similar to the moment vector but uses normalized central moments which take into account the frequency and distribution of nucleotides or amino acids from the origin of the sequence.

All of the Frequency Based and Natural Vector calculations are capable of being performed with different size k-nucleotides.  This value is adjustable by adjusting the Size of Kmers.

Graph Sequences:
	If the Moment Vector or Natural Vector buttons are selected two graphs are displayed: 
		1) A graphical representation of the DNA or amino acid sequence based on a vector system in which each of the four nucleotides corresponds to a different vector in the first or fourth quadrant of the Cartesian coordinate system.  Points in the graphical representation are obtained by summing the vectors representing nucleotides in the sequence.  The ordering of the nucleotide vectors is related to the GC content of genomes.  The amino acid system also uses two quadrants of the Cartesian coordinate system but orders its vectors according to the hydrophobicity of amino acids.
		2) A plot of a two-dimensional genome space using the components of the moment or natural vector (depending on which is selected) for each species.  The first component of the vector is plotted along the x-axis.  The y-axis plots another component of the vector based upon the selected dimension.  By default the selected dimension is 2.
	If one of the Frequency calculations is selected:
		A distance matrix is created for the selected species and the distances between the organisms are graphed in a 2D and 3D space.
		

Tree Type:
	GECKOV uses three distance matrix programs written by Joseph Felsenstein, University of Washington, as part of the Phylip package to draw a variety of phylogenetic trees based on various parameters.  A brief summary of the options available to the user in GECKOV are described below.  For a more complete description of the programs view the available documents in the README folder.

	Neighbor Joining: 
		Constructs a tree by successive clustering of lineages, setting branch lengths as the lineages join.  The tree is not rearranged thereafter.  The tree does not assume an evolutionary clock, so that it is effect an unrooted tree.  It should be somewhat similar to the tree obtained by Fitch.  The algorithm used is much faster than Fitch or Kitsch.  An outgroup should be specified.  The Jumble option does not allow multiple jumbles as there is no objective way of choosing which of the multiple results is best, there being no explicit criterion for optimality of the tree.  The Power option is not available for neighbor.  It is implicitly 0.
	UPGMA: 
		Unweighted Pair Group Method with Arithmetic Mean.  An option from within the Neighbor-Joining program, it is one of the simplest methods of tree construction.  It constructs a tree by successive clustering using an average-linkage method of clustering.  The tree topology will be the same as that of Kitsch.  The branch lengths of UPGMA will be the same as Kitsch with Power set at 0.  No outgroup option is available.
	Fitch: 
		Uses a weighted least squares method (Fitch-Margoliash) for clustering based on genetic distance.  The Minimum Evolution (ME) method is also available.  The ME method uses the Fitch-Margoliash criterion to fit branch lengths to each topology but then chooses topologies based on their total branch length (rather than the goodness of fit sum of squares).  There is no constraint on negative branch lengths in the ME method, so sometimes can give strange results.  This is the slowest of the algorithms as the time of the computation increases as n^4, where n is the number of species, compared to n^3 with the other algorithms.
	Kitsch: 
		Also uses the Fitch-Margoliash method, with the assumption that all tip species are contemporaneous and that there is an evolutionary clock.  The total length from the root of the tree to any species is the same.  No outgroup option is available as Kitsch estimates a rooted tree that cannot be rerooted.  Global rearrangements is set by default.

Global Rearrangements (G): 
	After the last species is added to the tree each possible group is removed and re-added improving the result since the position of every species is reconsidered.  Approximately triples the run-time.  Default in Kitsch.

Jumble (J):  
	Use a random number generator to choose the input order of the species.   Each different seed number leads to a different sequence of addition of species.  For a value of 10 (recommended), the program will try ten different orders of species in constructing the trees and the best tree from among all 10 runs will be chosen.

Power (P):
	The default value is 2 for the Fitch-Margoliash method.  This option is not available for Neighbor or UPGMA (both assume Power = 0).  The above methods assume that the variance of the measurement error is proportional to the P-th power of the expectation.  If you have reason to think that the measurement error of a distance is the same for small distances as it is for large, then you should set P=0 and use the least squares method, but if you have reason to think that the relative (percentage) error is more nearly constant than the absolute error, you should use P=2, the Fitch-Margoliash method. In between, P=1 would be appropriate if the sizes of the errors were proportional to the square roots of the expected distance.

Outgroup (O):
	This specifies which species is to have the root of the tree be on the line leading to it.  It is not available in UPGMA or Kitsch.  When it is used, the tree is still listed as being an unrooted tree, though the outgroup is connected to the bottommost node so that is easy to visually convert the tree into rooted form.

Retrieve Distances:
	Creates distance matrix and draws trees using the selected options for selected sequences from the box on the left.