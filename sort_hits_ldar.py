# Might be needed @cluster
#!/usr/local/bioinfo/src/python/Python-3.6.1/bin/python3.6
import os
import sys
import timeit

# Hard coded : set the correct path to edirect
edirect = '/home/scaonp01/Software/Edirect'

# Usage & warnings
usage = '\t --------\n' \
        '\t| usage  : python sort_hits_ldar.py f1\n' \
        '\t| input  : f1 = blastn.tsv (vs nt)\n' \
        '\t| stdout : % of reads <=> Salmonella, Cronobacter or others\n' \
        '\t --------'

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

##############
### Step 1 ###
##############
print('\n\tStep 1) Retrieve taxonomic infos for each staxids with efetch')
t0 = timeit.default_timer()

# For each line in TSV (f1), fill staxids_set
staxids_set = set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        columns = row.split('\t')
        # Sometimes you have more than one staxids for an entry
        staxids = columns[15].split(';')
        for i in staxids:
            staxids_set.add(i)

# Use staxids_set as query with efetch & store result
# Don't give more than let's say 500 entries at a time to avoid timeout
staxids_li = list(staxids_set)
staxids_sub_li = [staxids_li[x:x + 500]
                  for x in range(0, len(staxids_li), 500)]
efetch_li = []
for item in staxids_sub_li:
    staxids_input = ','.join(str(z) for z in item)
    # Details about "cmd" : https://www.biostars.org/p/163595/#271497
    cmd = (edirect + '/efetch -db taxonomy -id ' + staxids_input + ' -format '
           'xml | ' + edirect + '/xtract -pattern Taxon -sep \'@\' -element '
           'TaxId,ScientificName -division LineageEx -group Taxon -if Rank '
           '-equals superkingdom -or Rank -equals kingdom -or Rank -equals '
           'phylum -or Rank -equals class -or Rank -equals order -or Rank '
           '-equals family -or Rank -equals genus -sep \'@\' -element '
           'Rank,ScientificName')
    cmd_result = os.popen(cmd).read()
    cmd_result_split = cmd_result.split('\n')
    for i in cmd_result_split:
        efetch_li.append(i)

# Create a dict associating key=staxid with value=list=tax_infos
taxonomy_dic = {}
for line in efetch_li:
    field = line.split('\t')
    tax_ids = field[0].split('@')
    # Sometimes more than one staxid is associated to an entry
    # e.g. "170850@3666@Cucurbita hybrid cultivar"
    for i in tax_ids[:-1]:
        taxonomy_dic.setdefault(i, [None, None, None, None, None, None, None])
        for item in field:
            if 'superkingdom@' in item:
                taxonomy_dic[i][0] = item.split('@')[-1]
            elif 'kingdom@' in item:
                taxonomy_dic[i][1] = item.split('@')[-1]
            elif 'phylum@' in item:
                taxonomy_dic[i][2] = item.split('@')[-1]
            elif 'class@' in item:
                taxonomy_dic[i][3] = item.split('@')[-1]
            elif 'order@' in item:
                taxonomy_dic[i][4] = item.split('@')[-1]
            elif 'family@' in item:
                taxonomy_dic[i][5] = item.split('@')[-1]
            elif 'genus@' in item:
                taxonomy_dic[i][6] = item.split('@')[-1]

print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

##############
### Step 2 ###
##############
print('\tStep 2) Assign contigs best hit to Salmonella, Cronobacter or others')
t0 = timeit.default_timer()

# Assign contigs best hits to plant or non-plant based on taxonomy_dic infos
qseqid_set, salmo_hit_set, crono_hit_set, others_set = set(), set(), set(), set()
with open(sys.argv[1], 'r') as tsv:
    for row in tsv:
        columns = row.split('\t')
        qseqid, staxids = columns[0], columns[15].split(';')[0]
        # Check if we encounter qseqid for the first time <=> best hit
        if not qseqid in qseqid_set:
            if taxonomy_dic[staxids][6] == 'Salmonella':
                salmo_hit_set.add(qseqid)
            elif taxonomy_dic[staxids][6] == 'Cronobacter':
                crono_hit_set.add(qseqid)
            else:
                others_set.add(qseqid)
        qseqid_set.add(qseqid)

print('\t\t=> ' + str(round(timeit.default_timer() - t0, 3)) + ' seconds\n')

# Print results to stdout
print('\tSalmonella read count : ' + str(len(salmo_hit_set)))
print('\tCronobacter read count : ' + str(len(crono_hit_set)))
print('\tOthers read count : ' + str(len(others_set)) + '\n')
