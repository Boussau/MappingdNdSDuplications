# Introduction
mapDnDsSuplications is a program to map both substitution events and events of duplications and losses. in spite of the name, it can work not only with codons, but also with DNA or amino-acid alphabets. The syntax to the program follows the Bppsuite conventions. To run it, you need a reconciled gene tree in NHX format output by Phyldog, the corresponding species tree, and the alignment.

# Comiplation
A static compiled executable is provided for Linux systems. Otherwise I provide a basic compile command in the "compile" file.

# Output format
If site.and.branch=y, then a text file is produced, with lines such as: 
event(1,"test",L->M(2))
event(2,"test",M->N(3))
event(1,"test",duplication)
event(3,"test",loss)

This is telling us that there was a substitution from L to M on branch 1 of the species tree at site 2 of the alignment, and from M to N on branch 2 of the species tree at site 3. Similarly, the program found a gene duplication on branch 1, and a loss on branch 3. All of these events affected the gene family named "test".
To get the species tree node IDs and figure out where on the species tree those events occurred, a NHX species tree is output (speciesTreeName_withId), so that an integer is associated to each node of the species tree.

If site.and.branch=n, then Newick gene tree files are produced, with branch lengths corresponding to numbers of mapped events (e.g. one tre for dN events, one tree for dS events if map.type=DnDs).

# Examples to run the program: 

## Mapping all substitution events on protein sequences (quick and dirty, because the optimization is stopped very quickly, after 2 rounds or as soon as the loglikelihood improves by less than 10 points. It's probably better to use e.g. 1000 and 0.1 if you want to have robust substitution counts, but it does not matter if you're just interested in duplications and losses.). The option site.and.branch=y means that a table with events per site and per branch (including duplication/loss events) is output. Ancestral state reconstruction is used, meaning that 1 substitution is counted if there is a difference at a site between ML reconstructed parent and child sequences.
rm test_counts.txt ; ../mapDnDsDuplicationTree species.tree.file=test.Sptree input.tree.file=test.tree  input.sequence.file=test.fa  input.sequence.sites_to_use=all output.counts.file=test_counts.txt output.tree_with_id.file=test_id.nwk family.name=test map.type=All site.and.branch=y alphabet=Protein optimization.tolerance=10 optimization.max_number_f_eval=2 ancestral.state.reconstruction=y

## Mapping dN and dS substitutions on codon sequences, using the option site.and.branch=n, meaning that dN and dS numbers of substitutions are output as trees (one tree in which branch lengths are dN, one tree in which they are dS).
rm test_counts4.txt ; ../mapDnDsDuplicationTree species.tree.file=test4.Sptree input.tree.file=test4.tree  input.sequence.file=test4.fa  input.sequence.sites_to_use=all output.counts.file=test_counts4.txt output.tree_with_id.file=test_id4.nwk family.name=test4 map.type=DnDs site.and.branch=n 'alphabet=Codon(letter=DNA)' optimization.tolerance=10 optimization.max_number_f_eval=2 ancestral.state.reconstruction=y output.counts.tree.prefix=treeCounts4

## Same as above, but this time we do not use ancestral state reconstruction. This time we use proper substitution mapping, whcih means that at each site the expected number of substitutions of each type (here, dN and dS) are computed. Output as trees as above. This should be more robust than ancestral reconstruction based techniques, but is less useful when one wants to pinpoint precisely what transitions occurred at each site and on each branch.
rm test_counts4.txt ; ../mapDnDsDuplicationTree species.tree.file=test4.Sptree input.tree.file=test4.tree  input.sequence.file=test4.fa  input.sequence.sites_to_use=all output.counts.file=test_counts4.txt output.tree_with_id.file=test_id4.nwk family.name=test4 map.type=DnDs site.and.branch=n 'alphabet=Codon(letter=DNA)' optimization.tolerance=10 optimization.max_number_f_eval=2 ancestral.state.reconstruction=n output.counts.tree.prefix=treeCounts4


# Disclaimer and contact
This is a work in progress and most certainly contains some bugs, so please use with appropriate care and check that the results make sense. For bug reports and feature requests please contact me at boussau@gmail.com.
