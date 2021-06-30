#Divide and Cluster, currently just on the longest internal branch.

import sys, os
from ete3 import Tree
from Bio import SeqIO

def cut_on_longest_internal_branch(treeObj):
    branch_lengths = {}
    for node in treeObj.traverse("postorder"):
        branch_lengths[node] = node.dist
    sorted_branchlengths = sorted(branch_lengths.items(), key=lambda item: item[1], reverse=True)
    for node, branchlength in sorted_branchlengths:
        if node.is_leaf():
            continue
            #perhaps want to warn here that have some tips that are longer than the longest internal branch
        elif len(node) < 4:
            continue
        else:
            subtree1 = node.detach()
            subtree2 = treeObj
            return(subtree1, subtree2)
    return (treeObj, "No cut")

def create_subalignment(subtree, label): #given a subtree object, create a sequence set for those sequences, realign, retrim and run IQ-TREE bootstrapping
    clan_seqs = {}
    for tip in subtree:
        if tip.name in sequences:
            clan_seqs[tip.name] = str(sequences[tip.name].seq)
        else:
            print("Quitting: couldn't retrieve " + str(tip.name) + " from sequences file.")
            quit()
    clan_sequences_name = sys.argv[1] + "_" + label + "_sequences.fas"
    clan_alignment_name = sys.argv[1] + "_" + label + "_alignment.aln"
    subtree_name = sys.argv[1] + "_" + label + ".tre"
    outh = open(clan_sequences_name, "w")
    for seq in clan_seqs:
        outh.write(">" + seq + "\n" + clan_seqs[seq] + "\n")
    outh.close()
    os.system("mafft --auto " + clan_sequences_name + " > " + clan_alignment_name) #for speed in testing
    #os.system("mafft --localpair --maxiterate 1000 " + clan_sequences_name + " > " + clan_alignment_name) 
    os.system("divvier -partial -mincol 4 -divvygap " + clan_alignment_name)
    os.system("fasttree -lg -gamma " + clan_alignment_name + ".partial.fas > " + subtree_name) #for speed in testing, could add option to do this for the dividing at the start, too.



#main part    
sequences = SeqIO.index(sys.argv[1], "fasta")
tree = Tree(sys.argv[2])

(s1, s2) = cut_on_longest_internal_branch(tree)
if s2 == "No cut":
    print("Quitting: nowhere to divide the starting tree.")
    quit()

#now create sequence sets for the subtrees, align, and run trees
create_subalignment(s1, "clan1")
create_subalignment(s2, "clan2")



#trees_done = []
#trees_todo = []

#trees_todo.append(tree)
#while len(trees_todo) > 0:
#    (s1, s2) = bissect_tree(trees_todo.pop())
#    if s2 == 'Zilch': #no problem branches
#        trees_done.append(s1)
#    else:
#        trees_done.append(s1) #This relies on postorder and levelorder (gradually go deeper into tree), to prevent an infinite loop when cutting
#        trees_todo.append(s2)

#Write out the results
#counter = 0
#for tree in trees_done:
#    tree.write(outfile=sys.argv[2][:-4] + "_" + str(counter) + ".tre")
#    counter += 1
