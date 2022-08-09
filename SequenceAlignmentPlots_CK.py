## Code adapted from Adam Hockenberry's work from publication: D’Aquino,A.E., Azim,T., Aleksashin,N.A., Hockenberry,A.J., Krüger,A. and Jewett,M.C. (2020) Mutational characterization and mapping of the 70S ribosome active site. Nucleic Acids Res., 48, 2777–2789.
## available at https://github.com/adamhockenberry/23s-alignment-LTP


from Bio import Phylo
from Bio import SeqIO

import numpy as np
from scipy import stats
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42


## Names of libraries and unique sequences to find those locations (start of sequence matches start of library region, plus extra so that the search is specific)
name_list=['Helix 73 Segment 1, Residues 2046-2050','Helix 73 Segment 2, Residues 2618-2622','Helix 75 Segment 1, Residues 2088-2091','Helix 75 Segment 2, Residues 2228-2231','Helix 91, Residues 2523-2527' 'Helix 91, Residues 2536-2540','Helix 92: Residues 2457-2551','Helix 92: Residues 2557-2561']
target_seq_list=['GCGGCAAGACGGAAAGA','GCCGUGGGCGCUGGAGA','ACACUGAACAUUGAGC','GUGUCUGGUGGGUAGU','GGGGCUGAAGUAGGUCC','GUCCCAAGGGUAUGGCU','AUGGCUGUUCGCCAUUU','GCCAUUUAAAGUGGUAC']


for figtomake in range(0,len(name_list)):
    name=name_list[figtomake]
    target_seq=target_seq_list[figtomake]
    
    ## to define how long the library that needs to be plotted actually is. 
    if '73' in name:
        length=5
    if '75' in name:
        length=4
    if '91' in name:
        length=5
    if '92' in name:
        length=5
        
        
        
    print(name)
    
    ### read in tree data
    tree = Phylo.read('./LTPs123_LSU_tree.newick', 'newick')

    for i in tree.get_terminals()[:5]:
        print(i.name)

    ## read in sequences from fasta
    records = list(SeqIO.parse('./LTPs123_LSU_aligned.fasta', 'fasta'))

    ## find e coli
    for i in records:
        if i.description.find('Escherichia') != -1:
            print(i.description)


    sequence_dict = {}
    for record in records:
        sequence_dict[record.id] = record.seq

    print('There are {} sequences total and the length of the aligned E. coli sequence is {}'\
          .format(len(sequence_dict.keys()), len(sequence_dict['AJ278710'])))

    name_conversion_dict = {}
    for i in sequence_dict.keys():
    #     print(i)
        i_counts = 0
        for terminal in tree.get_terminals():
            if terminal.name.find(i) != -1:
                name_conversion_dict[terminal.name] = i
                i_counts += 1
        if i_counts != 1: ### Make sure that I find one and only one match for each record
            print(i, i_counts)



    mapping_dict = {}
    starting_point = 0
    nt_seq = []
    for i,nt in enumerate(list(str(sequence_dict['AJ278710']))):
        if nt in ['A', 'U', 'G', 'C']:
            nt_seq.append(nt)
            starting_point += 1
            mapping_dict[starting_point] = i

    ec_nt = ''.join(nt_seq)
    print('5\' region')
    print(ec_nt[:100])
    print('3\' region')
    print(ec_nt[-100:])


######################
    start=ec_nt.find(target_seq)

    start1=mapping_dict[start]+1
    end1=mapping_dict[start+len(target_seq)]

    for i in range(start, start+len(target_seq)+1):
        start1=mapping_dict[start]+1
        end1=mapping_dict[start+len(target_seq)]+1
        #print(mapping_dict[i])


    region_of_interest = list(range(start1, end1))


    start1_23Scode=start1-7015
    end1_23Scode=start1_23Scode+len(target_seq)
    start1_23Scode,end1_23Scode

    residue_numbering=[]
    for i in range(start1_23Scode,end1_23Scode):
        residue_numbering.append(i)
    residue_numbering

    region_name=name#+': residue '+ str(start1_23Scode)+ ' to '+ str(end1_23Scode)
    region_name
    
    ## in case there are gaps in the alignment and we don't want to plot the gaps. If you do want to include the gaps, 

    default='NoGaps'
    deletionloc=default
    residue_ids=([sequence_dict['AJ278710'][i] for i in region_of_interest])
    if '-' in residue_ids:
        deletionloc=residue_ids.index('-')

    matrix = []

    for terminal in tree.get_terminals():
        i = name_conversion_dict[terminal.name]
        j = sequence_dict[i]
        tempy = []
        if len(j) > 1600:
            for nt_pos in range(0,len(region_of_interest)):#region_of_interest:
                if deletionloc!=default:
                    if nt_pos!=deletionloc:
                        nt = j[region_of_interest[nt_pos]]
                        #if nt !='-': deletionloc
                        if nt == 'U':
                            tempy.append(1)
                        elif nt == 'C':
                            tempy.append(0.5)
                        elif nt == 'A':
                            tempy.append(0)
                        elif nt == 'G':
                            tempy.append(-0.5)
                        else:
                            tempy.append(-1)
            matrix.append(tempy)

    print('Double checking that my matrix has {} rows and {} columns'.format(len(matrix), len(matrix[0])))

    ### Plottign the color map.
    region_of_interest_clean=[]
    for position in range(0,len(residue_ids)):
        if residue_ids[position]!='-':
            region_of_interest_clean.append(residue_ids[position])

            
    matrix_edit=[i[0:length] for i in matrix] ## this is to take only our library regions instead of the whole search sequence

    fig, ax = plt.subplots(figsize=(3+(9*(len(matrix_edit[0])/22)),12))

    # define the colormap
    cmap = plt.cm.jet
    cmaplist = ['white', 'lightseagreen', 'rebeccapurple', 'darkorange', 'royalblue']

    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(-1.2, 1.2, 6)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    cax = ax.pcolormesh(matrix_edit, norm=norm, cmap=cmap)
    cax.set_edgecolor('face')
    cbar = fig.colorbar(cax, ticks=[-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(['Missing', 'G', 'A', 'C', 'U'], fontsize=24)
    ax.set_xticks([i+0.5 for i in range(len(matrix_edit[0]))])
    ax.set_xticklabels(region_of_interest_clean)
    ax.xaxis.tick_top()
    ax.set_xlabel('NT position', fontsize=24)
    ax.set_ylabel('Species', fontsize=24)
    ax.tick_params(labelsize=20)
    ax.set_yticklabels('')

    for line in ax.yaxis.get_ticklines():
        line.set_markersize(0)
        line.set_markeredgewidth(0)
    plt.tight_layout()
    plt.title(region_name)
    plt.savefig(name+'.pdf', bbox_inches='tight')


     #NOW MAKING THE SEQ CONSERVATION BARPLOTS!

    freqdict_total=pd.DataFrame({})
    freqdict_total['Residue Identity']=['G','A','C','U']
    #freqdict_total.set_index('Residue Identity')


    from collections import Counter

    percent_conserved=[]
    freqs_log=[]
    counterlog=[]
    for pos, i in enumerate(np.array(matrix).T):
       # if pos !=deletionloc:
        if pos<len(region_of_interest_clean):
            save=(pos,i)

            vallist=(save[1])
            emptycount=[]
            Gcount=[]
            Acount=[]
            Ccount=[]
            Ucount=[]
            for i in range(0,len(save[1])):
                if vallist[i]==-1:
                    emptycount.append(vallist[i])
                if vallist[i]==-0.5:
                    Gcount.append(vallist[i])
                if vallist[i]==0:
                    Acount.append(vallist[i])
                if vallist[i]==0.5:
                    Ccount.append(vallist[i])
                if vallist[i]==1:
                    Ucount.append(vallist[i])
                freq_dict={'G':(len(Gcount)/len(save[1])),'A':(len(Acount)/len(save[1])),'C':len(Ccount)/len(save[1]),'U':len(Ucount)/len(save[1])}
            freqdict_total[pos]=list(freq_dict.values())
            
            if pos<len(region_of_interest_clean):
                percent_conserved.append(freq_dict[region_of_interest_clean[pos]])

    ## to save conservation data to excel
    pd.DataFrame({'Residue':region_of_interest_clean[0:length], '% Conservation':percent_conserved[0:length]}).to_excel(name+'_conservation.xlsx')


    freqdict_total=freqdict_total.set_index('Residue Identity')

    cmaplist = ['white', 'lightseagreen', 'rebeccapurple', 'darkorange', 'royalblue']
    cmaplist_noMissing=cmaplist[1:]

    transposed=freqdict_total.transpose()
    transposed['ResList']=region_of_interest_clean
    new=transposed.set_index('ResList')
    new=new[0:length]

    import seaborn as sns
    new.plot(kind='bar',color=cmaplist_noMissing,width = 0.7,figsize=(len(region_of_interest_clean),5))
    #plt.figure(figsize=(0.5*len(residue_ids),5))
    plt.xticks(rotation=0)
    plt.grid(False)
    plt.xlabel('Residue Identities')
    plt.ylabel('% Conservation')
    plt.title(name)
    plt.legend(bbox_to_anchor=(1,1))
    plt.savefig(name+'_barplot.pdf')