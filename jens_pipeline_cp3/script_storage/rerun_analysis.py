from __future__ import division
from __future__ import print_function
from operator import itemgetter
import argparse
import os
import sys
import shutil
import subprocess
import itertools
import time

##Version 0.9.8.3 (Increased low complex search span, fixed papa number in candidate list)

# Description:
# This scripts takes the output summary files from all different prediction softwares
# and compiles them into a single result. This script is intended to rerun the standard 
# pipeline final analysis with varying settings.
# The first sweep selects candidates based on co-prediction of the intrinsic disorder softwares, 
# the plaac score and papa score. The second sweep is designed to incorporate low complexity
# in the analysis. This is done in this fashion because the low complexity analysis
# bottlenecks the pipeline quite severely, hence it is only run on the first round 
# candidates. Candidates are thus proteins whose vitarbi parse (with sufficient LLR-score)
# are found to overlap or contain segments of intrinsic disorder as judged by majority rule
# (at least 3 out of 4 potential softwares) as well as low complexity, whilst also 
# having a high enough papa.py score.
# Candidates are thus proteins that contain multiple qualities that each contribute to 
# aggregation and/or amyloidal structuring.


# Default values:
# args.filter_disorder : 9 (removes below and including)
# args.filter_low_complexity : 4 (removes below and including)
# args.disorder: 3 (max 4)
# args.low_complexity : False
# args.plaac : >=10 (somewhat lenient)
# papa_score_cutoff : >= 0.05 (program default)

# Parser for working flags to be used in command line.
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",type=str,help="Input folder name")
parser.add_argument("-o","--output",type=str,help="Output folder_name")
parser.add_argument("-t","--tag",type=str,help="tag name used for analysis")
parser.add_argument("-pl","--plaac",type=int,help="plaac score cut-off value. Default 10",default=10)
parser.add_argument("-d","--disorder",type=int,\
                    help="Required number of disorder programs to co-predict. Default 4",default=3)
parser.add_argument("-fid","--filter_disorder",type=int,\
                    help="Cut-off for filtering disorder. Default 9 and below",default=9)
parser.add_argument("-flc","--filter_low_complexity",type=int,\
                    help="Cut-off for filtering low-complexity. Default 4 and below",default=4)
parser.add_argument("-lc","--low_complexity",type=str,\
                    help="Include low-complexity in prion analysis. Default False",\
                    default="False")
parser.add_argument("-pa","--papa",type=str,\
                    help="Include papa.py in analysis. Default True",\
                    default="True")
args=parser.parse_args()

yesses=["t","true","yes","y"]
nosses=["f","false","no","n"]

# Checks that all required flags are present
if not args.input or not args.tag or not args.output:
    print("INPUT ERROR: Required flag missing: -i,-o,-t",file=sys.stdout)
    sys.exit()

# Checks that the string input for toggling low complexity is one of the accepted answers.
if args.low_complexity.lower() not in yesses and args.low_complexity.lower() not in nosses:
    print("Error: Incorrect input for -lc. Accepted answers include [yes,no,true,false]",file=sys.stdout)
    sys.exit()
if args.low_complexity.lower() in yesses:
    args.low_complexity = True
else:
    args.low_complexity = False

# Checks that the string input for toggling papa.py is one of the accepted answers.
if args.papa.lower() not in yesses and args.papa.lower() not in nosses:
    print("Error: Incorrect input for -pa. Accepted answers include [yes,no,true,false]",file=sys.stdout)
    sys.exit()
if args.papa.lower() in nosses:
    args.papa = False
else:
    args.papa = True

base_path = os.getcwd()


#Ensures the correct amount of slashes for file path
if args.input[-1] == "/":
    path_to_input = "../" + str(args.input)
else:
    path_to_input = "../" + str(args.input) + "/"

#Checks if output directory name already exists
if os.path.exists(args.output):
    print("------------------------------------------------------------",file=sys.stdout)
    print("INPUT ERROR:",file=sys.stdout)
    print("Output directory name already used, please select a new one:",file=sys.stdout)
    sys.exit()

# Function reading intrinsic disorder span files, writing 'protein name : span ranges' to dictionary
def disorder(file_lines, filter_length):
    regions = []
    name = ""
    span_dict = {}

    for x in file_lines:
        # Continue if line is empty
        if x[0] == "\n":
            continue

        # Name is set according to identifier '>'
        if x.startswith(">"):
            x = x.rstrip("\n")
            name = x[1:]

        # If end of span for protein, append placeholder regions to full dictionary
        if x.startswith("#"):
            regions = [j for i in regions for j in i]
            span_dict[name] = regions
            name = ""
            regions = []

        # If line is not identifier or end, append number to placeholder 'regions'
        if x[0] != ">" and x[0] != "#" and x[0] != "\n":
            x = x.rstrip("\n")
            x = x.split('-')
            y = int(x[1])
            z = int(x[0])
            x = y + 1 - z

            # If length of span below filter length, simply ignore the span and continue
            if filter_length > 0:
                if x <= filter_length:
                    continue
                else:
                    regions.append(range(z, y + 1))
            else:
                regions.append(range(z, y + 1))
        else:
            continue
    return span_dict

# Read all result files from pipeline
with open("iup_short_spans", "r") as a, open("vsl2_spans", "r") as b, \
        open("espritz_x_spans", "r") as c, open("espritz_n_spans", "r") as d, \
        open("papa_summary", "r") as e, open("plaac_output", "r") as f:
    iup_s = a.readlines()
    vsl2 = b.readlines()
    esp_x = c.readlines()
    esp_n = d.readlines()
    papa = e.readlines()
    plaac = f.readlines()

# Send the list of lines to function 'disorder', write dictionary to variable.
iup_s_dict = disorder(iup_s, args.filter_disorder)
vsl2_dict = disorder(vsl2, args.filter_disorder)
esp_x_dict = disorder(esp_x, args.filter_disorder)
esp_n_dict = disorder(esp_n, args.filter_disorder)

# Read lines from papa-file and create dictionary, format 'Name:[score,position]'
# Ignores lines with string "protein length too short..."
papa_dict = {}
for line in papa:
    if "p" in str(line.split("\t")):
        continue
    line = line.strip()
    line = line.split("\t")

    papa_dict[str(line[0])] = [float(line[1]), int(line[2])]

# Read lines from plaac-file and create dictionary, format 'Name:[score,middle position, vitarbi parse string,[middle position +- 25]'
plaac_dict = {}
for line in plaac[1:]:
    line = line.strip()
    line = line.split("\t")
    middle_llr = int((int(line[6]) + int(line[7])) / 2)
    lower_llr = middle_llr-25
    upper_llr = middle_llr+26
    plaac_dict[line[0]] = [float(line[11]), middle_llr, str(line[22]),range(lower_llr,upper_llr)]


# Create list of candidates according to number of agreeing intrinsic disorder
# programs, and a sufficient plaac and papa-score.
candidates = []
for key, value in plaac_dict.iteritems():
    score_id = 0
    if value[0] >= args.plaac and value[0] != "NaN" and args.papa==False:
        if value[1] in iup_s_dict[key]:
            score_id += 1
        if value[1] in vsl2_dict[key]:
            score_id += 1
        if value[1] in esp_x_dict[key]:
            score_id += 1
        if value[1] in esp_n_dict[key]:
            score_id += 1

    elif value[0] >= args.plaac and value[0] != "NaN" and args.papa==True:
        if papa_dict[key][0] > 0.05:
            if value[1] in iup_s_dict[key]:
                score_id += 1
            if value[1] in vsl2_dict[key]:
                score_id += 1
            if value[1] in esp_x_dict[key]:
                score_id += 1
            if value[1] in esp_n_dict[key]:
                score_id += 1
        else:
            continue
    else:
        continue

    if score_id >= args.disorder:
        candidates.append(key)
    else:
        continue

# Checks if the candidate list is empty. Aborts remaining scan if empty.
if len(candidates) == 0:
    print("No primary candidates found, aborting scan...",file=sys.stdout)
    sys.exit()
# If there are any candidate proteins, and low complexity is toggled off,
# create the file "candidate_prions..." and append all candidates with their
# plaac score to this file. This file is later read by other sections of the script.
# Also sorts the output according to decreasing plaac score.
if len(candidates) > 0 and args.low_complexity == False:
    os.makedirs(args.output)
    os.chdir(args.output)
    temp_sort_list=[]
    g = open("candidate_prions_first_sweep", "w+")
    for protein in candidates:
        temp_sort_list.append([protein,plaac_dict[protein][0],papa_dict[protein][0]])
    temp_sort_list.sort(key=itemgetter(1), reverse=True)
    print("Systematic Name:\tPlaac LLR-score:\tPapa score:", file=g)
    for id in temp_sort_list:
        print(id[0]+"\t"+str(id[1])+"\t"+str(id[2]), file=g)
    g.close()


# Same as above, but with low complexity toggled on. With the created candidate list,
# send the sequences of those candidate proteins to GBA for low complexity analysis.
# Then, retrieve the GBA summary file, and delete the sequences in the GBA folder.
if len(candidates) > 0 and args.low_complexity == True:
    os.makedirs(args.output)
    os.chdir(args.output)
    temp_sort_list=[]
    g = open("candidate_prions_first_sweep", "w+")
    for protein in candidates:
        temp_sort_list.append([protein,plaac_dict[protein][0],papa_dict[protein][0]])
    temp_sort_list.sort(key=itemgetter(1),reverse=True)
    print("Systematic Name:\tPlaac LLR-score:\tPapa score:",file=g)
    for id in temp_sort_list:
        print(id[0]+"\t"+str(id[1])+"\t"+str(id[2]), file=g)
    g.close()
    

    folder_gba = '../../programs/GBA/'

    with open("candidate_prions_first_sweep") as h:
        names = h.readlines()

    os.makedirs(folder_gba + args.output)

    for x in names[1:]:
        x = x.split(' ')
        x = x[0].strip()
        filename = x + ".fasta"
        shutil.copy("../" + path_to_input + filename, folder_gba + args.output + "/")

    subprocess.call([folder_gba + './pipeline_gba.sh', args.output],cwd=folder_gba)

    shutil.move(folder_gba + args.output + "_lcr_spans", os.getcwd())

# Read the GBA summary file containing low complexity spans into a dictionary.
# Pretty much the same format as for the intrinsic disorder dictionaries.
if args.low_complexity == True:
    with open(args.output+"_lcr_spans") as lcr_spans:
        lcr_lines = lcr_spans.readlines()

    gba_dict = disorder(lcr_lines, args.filter_low_complexity)

# Function for creating the final list of ranges of low complexity for each protein.
# Filters length that are below the set limit.
def lcr_extract(lcr_spans, protein, filter_len, prot_len):
    logic_gate = False
    span_holder = []
    temp_span_range = []
    binaries = ""
    for row in lcr_spans:
        row = row.strip("\n")
        if logic_gate == True and row[0] != "#":
            split_row = row.split("-")
            residue_diff = int(split_row[1]) - int(split_row[0]) + 1
            temp_span_range.extend(range(int(split_row[0]), int(split_row[1]) + 1))
            if int(residue_diff) <= filter_len:
                continue
            else:
                span_holder.append(str(row))
        if logic_gate == True and row[0] == "#":
            logic_gate = False
        if row.strip(">") == protein:
            logic_gate = True
    for residue in range(1, prot_len + 1):
        if residue in temp_span_range:
            binaries += "X"
        else:
            binaries += "-"

    summ_list = [span_holder, binaries]
    return summ_list

# Function for retrieving protein name, sequence length and creating a binary map
# for each protein, which is later written to the final summary result files.
def read_write(protein_file, filter_len):
    split_spans = []
    seq_length = 0
    binaries = ""
    temp_span_ranges = []
    for row in protein_file[1:]:
        if row[0] == "#":
            if int(row.split(":")[1]) <= filter_len:
                continue
            else:
                row = row.split(':')[0]
                row = row.strip('#')
                row = row.split('-')
                split_spans.append(int(row[0]))
                split_spans.append(int(row[1]))

        if row[0:8] == "Sequence":
            seq_length = int(row.split(': ')[1])

    residue_counter = 0
    for span_boundary in split_spans:
        if residue_counter % 2 == 0:
            lower = int(split_spans[residue_counter])
            upper = int(split_spans[residue_counter + 1])
            temp_span_ranges.extend(range(lower, upper))
        residue_counter += 1

    for residues in range(1, seq_length + 1):
        if residues in temp_span_ranges:
            binaries += "D"
        else:
            binaries += "-"
    span_couples = []
    split_spans_dummy = split_spans
    for x in split_spans_dummy:
        try:
            span_couples.append(str(split_spans_dummy[0]) + '-' + str(split_spans_dummy[1]))
            split_spans_dummy = split_spans_dummy[2:]
        except IndexError:
            continue

    summ_list = [span_couples, seq_length, binaries]
    return summ_list

# Makes a new folder to contain the first candidates summary files.
os.makedirs("first_sweep_summaries")

# For each candidate, and if low complexity is toggled on, create summary files
# that are written in the same manner as those without low complexity, only with
# the low complexity analysis included.
# Creates a second sweep of candidates according to low complexity,
# and moves result files to a secondary sweep result folder, as well as a
# new candidate list.
if args.low_complexity == True:
    write_dictionary = {}
    for protein in candidates:
        name = str(protein)
        protein = protein.split(' ')
        i = open("../iupred_short_files/iup_" + args.tag + "." + protein[0])
        i_lines = i.readlines()
        j = open("../vsl2_files/vsl2_" + args.tag + "." + protein[0])
        j_lines = j.readlines()
        k = open("../espritz_n_files/espr_" + args.tag + "_n." + protein[0])
        k_lines = k.readlines()
        l = open("../espritz_x_files/espr_" + args.tag + "_x." + protein[0])
        l_lines = l.readlines()
        lcr_span_file = open(args.output+"_lcr_spans")
        lcr_lines = lcr_span_file.readlines()

        iup_cand_write = read_write(i_lines, args.filter_disorder)
        vsl2_cand_write = read_write(j_lines, args.filter_disorder)
        espr_n_cand_write = read_write(k_lines, args.filter_disorder)
        espr_x_cand_write = read_write(l_lines, args.filter_disorder)
        prot_len = iup_cand_write[1]
        gba_cand_write = lcr_extract(lcr_lines, name, args.filter_low_complexity, prot_len)

        # Makes a master span list that itertools.izip_longest uses.
        master_span_list = [iup_cand_write[0], vsl2_cand_write[0], espr_n_cand_write[0], \
                            espr_x_cand_write[0], gba_cand_write[0]]

        # Prints all info into summary files from all the various dictionaries.
        o = open("summary." + str(protein[0]), "w+")
        print("Protein name:", name, file=o)
        print("Sequence length:", iup_cand_write[1], file=o)
        print("Plaac-score:", plaac_dict[name][0], file=o)
        print("Vitarbi parse:", plaac_dict[name][2], file=o)
        print("Vitarbi centre residue:", plaac_dict[name][1], file=o)
        print("Papa-score:", papa_dict[name][0], file=o)
        print("Papa-score position:", papa_dict[name][1], file=o)
        print("================================================================================", file=o)
        print("Filtering value ID:", args.filter_disorder, file=o)
        print("Filtering value LC:", args.filter_low_complexity, "\n", file=o)
        print("IUPred\tVSL2b\tEspritz N\tEspritz X\tGBA", file=o)
        print("________________________________________________________________________________", file=o)
        for x in itertools.izip_longest(*master_span_list, fillvalue=""):
            print('\t'.join([str(e) for e in x]), file=o)
        print("================================================================================", file=o)
        print("IUPred binary:", iup_cand_write[2], file=o)
        print("VSL2b binary:", vsl2_cand_write[2], file=o)
        print("Espritz N binary:", espr_n_cand_write[2], file=o)
        print("Espritz X binary:", espr_x_cand_write[2], file=o)
        print("GBA binary:", gba_cand_write[1], file=o)
        o.close()

    os.makedirs("second_sweep_summaries")

    # Reads all protein names from the first sweep candidate file, and sends those 
    # proteins for low complexity analysis. Proteins whose disorder and viterbi parse
    # agree with its low complexity are written to a second sweep file, and moved
    # to a second sweep summary directory.
    with open("candidate_prions_first_sweep") as first_sweep:
        second_sweep_file = open("candidate_prions_second_sweep", "w+")
        print("Systematic Name:\tPlaac LLR-score:\tPapa score:",file=second_sweep_file)
        first_sweep_rows = first_sweep.readlines()
        temp_sort_list=[]
        for row in first_sweep_rows[1:]:
            row=row.split('\t')
            prot_name = row[0].strip()
            syst_name = prot_name.split(' ')[0]
            shutil.move(os.getcwd() + "/summary." + syst_name, os.getcwd() + "/first_sweep_summaries")
            for residue in plaac_dict[prot_name][3]:
                if residue in gba_dict[prot_name]:
                    temp_sort_list.append([prot_name,plaac_dict[prot_name][0],papa_dict[prot_name][0]])
                    shutil.copy(os.getcwd() + "/first_sweep_summaries/summary." + syst_name, \
                            os.getcwd() + "/second_sweep_summaries")
                    break
                else:
                    continue
        temp_sort_list.sort(key=itemgetter(1),reverse=True)
        for id in temp_sort_list:
            print(id[0]+"\t"+str(id[1])+"\t"+str(id[2]),file=second_sweep_file)

# Same as above, except without low complexity. Does not create secondary sweep files.
if args.low_complexity == False:
    write_dictionary = {}
    for protein in candidates:
        name = str(protein)
        protein = protein.split(' ')
        i = open("../iupred_short_files/iup_" + args.tag + "." + protein[0])
        i_lines = i.readlines()
        j = open("../vsl2_files/vsl2_" + args.tag + "." + protein[0])
        j_lines = j.readlines()
        k = open("../espritz_n_files/espr_" + args.tag + "_n." + protein[0])
        k_lines = k.readlines()
        l = open("../espritz_x_files/espr_" + args.tag + "_x." + protein[0])
        l_lines = l.readlines()

        iup_cand_write = read_write(i_lines, args.filter_disorder)
        vsl2_cand_write = read_write(j_lines, args.filter_disorder)
        espr_n_cand_write = read_write(k_lines, args.filter_disorder)
        espr_x_cand_write = read_write(l_lines, args.filter_disorder)
        prot_len = iup_cand_write[1]


        # Makes a master span list that itertools.izip_longest uses.
        master_span_list = [iup_cand_write[0], vsl2_cand_write[0], espr_n_cand_write[0], \
                            espr_x_cand_write[0]]

        # Prints all info into summary files from all the various dictionaries.
        o = open("summary." + str(protein[0]), "w+")
        print("Protein name:", name, file=o)
        print("Sequence length:", iup_cand_write[1], file=o)
        print("Plaac-score:", plaac_dict[name][0], file=o)
        print("Vitarbi parse:", plaac_dict[name][2], file=o)
        print("Vitarbi centre residue:", plaac_dict[name][1], file=o)
        print("Papa-score:", papa_dict[name][0], file=o)
        print("Papa-score position:", papa_dict[name][1], file=o)
        print("=======================================================", file=o)
        print("Filtering value ID:", args.filter_disorder, file=o)
        print("IUPred\tVSL2b\tEspritz N\tEspritz X", file=o)
        print("_______________________________________________________", file=o)
        for x in itertools.izip_longest(*master_span_list, fillvalue=""):
            print('\t'.join([str(e) for e in x]), file=o)
        print("=======================================================", file=o)
        print("IUPred binary:", iup_cand_write[2], file=o)
        print("VSL2b binary:", vsl2_cand_write[2], file=o)
        print("Espritz N binary:", espr_n_cand_write[2], file=o)
        print("Espritz X binary:", espr_x_cand_write[2], file=o)

        o.close()

    with open("candidate_prions_first_sweep") as first_sweep:
        first_sweep_rows = first_sweep.readlines()
        for row in first_sweep_rows[1:]:
            prot_name = row.strip()
            syst_name = row.split(' ')[0]
            shutil.move(os.getcwd() + "/summary." + syst_name, \
                        os.getcwd() + "/first_sweep_summaries")

# Finally writes a file containing the settings used for that round of analysis, moved
# to the output folder specified from the start.
sum_file=open("settings_%s" %args.output,"w+")
print("Input folder:\t%s" %args.input,file=sum_file)
print("Output folder:\t%s" %args.output,file=sum_file)
print("Tag name:\t%s" %args.tag,file=sum_file)
print("Plaac cut-off:\t%s" %args.plaac,file=sum_file)
print("Nr disorder programs:\t%s" %args.disorder,file=sum_file)
print("Disorder filtering value:\t%s" %args.filter_disorder,file=sum_file)
print("Low complexity included:\t%s" %args.low_complexity,file=sum_file)
print("Low-complexity filtering value:\t%s" %args.filter_low_complexity,file=sum_file)
print("Papa.py included:\t%s" %args.papa,file=sum_file)
print("Date & time: " + time.strftime("%c"),file=sum_file)
sum_file.close()

