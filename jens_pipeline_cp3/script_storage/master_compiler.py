from __future__ import division
from __future__ import print_function
from operator import itemgetter
import os
import sys
import shutil
import subprocess
import itertools

# Master compiling script version 0.9.9 (include LC as default)
# Part of the default pipeline, initiated by the start_script.sh,
# compiles results of all predictors and creates lists of candidates
# from the results.

#struct: /file/to/path
input_folder=sys.argv[1]
tag_name=sys.argv[2]
output_folder=sys.argv[3]
base_path=sys.argv[4]

path_to_input=base_path+'/'+input_folder+'/'

os.chdir(base_path+'/'+output_folder+'/')

idp_filter_length=9
gba_filter_length=4



##Function reading span files, writing 'protein name : span ranges' to dictionary
def disorder(file_lines,filter_length):
        regions=[]
        name=""
        span_dict={}

        for x in file_lines:
            #Continue if line is empty
            if x[0]=="\n":
                continue
            
            #Name is set according to identifier '>'
            if x.startswith(">"):   
                x=x.rstrip("\n")
                name=x[1:]
                
            #If end of span for protein, append placeholder regions to full dictionary    
            if x.startswith("#"):
                regions=[j for i in regions for j in i]
                span_dict[name]=regions
                name=""
                regions=[]
                
            #If line is not identifier or end, append number to placeholder 'regions'   
            if x[0]!= ">" and x[0] != "#" and x[0] != "\n":
                x=x.rstrip("\n")
                x=x.split('-')
                y=int(x[1])
                z=int(x[0])
                x=y+1-z
                
                #If length of span below filter length, simply ignore the span and continue
                if filter_length > 0: 
                    if x <= filter_length:
                        continue
                    else:
                        regions.append(range(z,y+1))
                else:
                    regions.append(range(z,y+1))
            else:
                continue
        return span_dict

#Read all result files from pipeline
with open("iup_short_spans","r") as a, open("vsl2_spans","r") as b,\
     open("espritz_x_spans","r") as c, open("espritz_n_spans","r") as d,\
     open("papa_summary","r") as e, open("plaac_output","r") as f:
          iup_s=a.readlines()
          vsl2=b.readlines()
          esp_x=c.readlines()
          esp_n=d.readlines()
          papa=e.readlines()
          plaac=f.readlines()


#Send the list of lines to function 'disorder', write dictionary to variable.
iup_s_dict=disorder(iup_s,idp_filter_length)
vsl2_dict=disorder(vsl2,idp_filter_length)
esp_x_dict=disorder(esp_x,idp_filter_length)
esp_n_dict=disorder(esp_n,idp_filter_length)


#Read lines from papa-file and create dictionary, format 'Name:[score,position]'
papa_dict={}
for line in papa:
        if "p" in str(line.split("\t")):
                continue
        line=line.strip()
        line=line.split("\t")
        
        papa_dict[str(line[0])]=[float(line[1]),int(line[2])]

#Read lines from plaac-file and create dictionary, format 'Name:[score,middle position, vitarbi parse string]'
plaac_dict={}
for line in plaac[1:]:
        line=line.strip()
        line=line.split("\t")
        middle_llr=int((int(line[6])+int(line[7]))/2)
        plaac_dict[line[0]]=[float(line[5]),middle_llr,str(line[22])]

##Create dictionary of proteins according to number of intrinsic
##disorder programs co-predicting and coinciding with plaac prediction.
##Disregards proteins with plaac score below 10.
candidates = []
for key, value in plaac_dict.iteritems():
    score_id = 0
    if value[0] >= 10 and value[0] != "NaN" and papa_dict[key][1] > 0.05:
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

    if score_id >= 4:
        candidates.append(key)
    else:
        continue

#Checks if there are any candidates, and sends them to low complexity analysis.
#Sends result file from gba back to master summary directory.
if len(candidates) > 0:
    print(candidates)
    temp_sort_list=[]
    g = open("candidate_prions_first_sweep_def","w+")
    for protein in candidates:
        temp_sort_list.append([protein,plaac_dict[protein][0],papa_dict[protein][0]])
    temp_sort_list.sort(key=itemgetter(1),reverse=True)
    print("Systematic Name:\tPlaac LLR-score:\tPapa score:",file=g)
    for id in temp_sort_list:
        print(id[0]+"\t"+str(id[1])+"\t"+str(id[2]), file=g)
    g.close()

else:
    print("nothing!")
    
#
#
#    folder_gba= '../programs/GBA/'
#
#    with open("candidate_prions_first_sweep_def") as h:
#        names=h.readlines()
#
#    os.makedirs(folder_gba + sys.argv[3])
#    for x in names[1:]:
#        x=x.split(' ')
#        x=x[0].strip()
#        filename=x + ".fasta"
#        shutil.copy(path_to_input + filename, folder_gba + sys.argv[3] + "/")
#    subprocess.call([base_path + "/programs/GBA/./pipeline_gba.sh", str(sys.argv[3])],cwd=folder_gba)
#    shutil.move(base_path + "/programs/GBA/" + sys.argv[3] + "_lcr_spans", os.getcwd())

##Opens newly created LC span file and reads into a LC dictionary.
#with open(sys.argv[3] + "_lcr_spans") as lcr_spans:
#    lcr_lines=lcr_spans.readlines()
        
#gba_dict=disorder(lcr_lines,gba_filter_length)

##Function similar to the one below except specialized for low complexity.
##(Reads from different sources and are hence different, but same output)
def lcr_extract(lcr_spans,protein,filter_len,prot_len):
        logic_gate=False
        span_holder=[]
        temp_span_range=[]
        binaries=""
        for row in lcr_spans:
            row=row.strip("\n")
            if logic_gate==True and row[0]!="#":
                    split_row=row.split("-")
                    residue_diff = int(split_row[1])-int(split_row[0])+1
                    temp_span_range.extend(range(int(split_row[0]),int(split_row[1])+1))
                    if int(residue_diff) <= filter_len:
                            continue
                    else:
                            span_holder.append(str(row))
            if logic_gate==True and row[0]=="#":
                    logic_gate=False
            if row.strip(">")==protein:
                    logic_gate=True
        for residue in range(1,prot_len+1):
                if residue in temp_span_range:
                        binaries+="X"
                else:
                        binaries+="-"
        
        summ_list=[span_holder,binaries]
        return summ_list
                
##Function for extracting info being written to summary files for individual proteins.
##Extracts name, length, binary maps, disorder spans
def read_write(protein_file,name,filter_len):
        split_spans=[]
        seq_length=0
        binaries=""
        temp_span_ranges=[]
        for row in protein_file[1:]:
                if row[0]=="#":
                        if int(row.split(":")[1])<=filter_len:
                                continue
                        else:
                                row=row.split(':')[0]
                                row=row.strip('#')
                                row=row.split('-')
                                split_spans.append(int(row[0]))
                                split_spans.append(int(row[1]))

                if row[0:8]=="Sequence":
                        seq_length=int(row.split(': ')[1])
        
        residue_counter=0
        for span_boundary in split_spans:
                if residue_counter % 2==0:
                        lower=int(split_spans[residue_counter])
                        upper=int(split_spans[residue_counter+1])
                        temp_span_ranges.extend(range(lower,upper))
                residue_counter+=1
                                              
        for residues in range(1,seq_length+1):
                if residues in temp_span_ranges:
                        binaries+="D"
                else:
                        binaries+="-"
        span_couples=[]
        split_spans_dummy=split_spans
        for x in split_spans_dummy:
                try:
                        span_couples.append(str(split_spans_dummy[0])+'-'+str(split_spans_dummy[1]))
                        split_spans_dummy=split_spans_dummy[2:]
                except IndexError:
                        continue
                
        summ_list=[span_couples,seq_length,binaries]
        return summ_list

##Reads info from files according to candidate prion list created earlier.
##Opens files of each protein from each program and sends file rows to read_write function.
write_dictionary={}
for protein in candidates:
        name=str(protein)
        protein=protein.split(' ')
        i = open("iupred_short_files/iup_"+tag_name+"."+protein[0])
        i_lines=i.readlines()
        j = open("vsl2_files/vsl2_"+tag_name+"."+protein[0])
        j_lines=j.readlines()
        k = open("espritz_n_files/espr_"+tag_name+"_n."+protein[0])
        k_lines=k.readlines()
        l = open("espritz_x_files/espr_"+tag_name+"_x."+protein[0])
        l_lines=l.readlines()
        #lcr_span_file= open(sys.argv[3] + "_lcr_spans")
        #lcr_lines=lcr_span_file.readlines()

        iup_cand_write=read_write(i_lines,name,idp_filter_length)
        vsl2_cand_write=read_write(j_lines,name,idp_filter_length)
        espr_n_cand_write=read_write(k_lines,name,idp_filter_length)
        espr_x_cand_write=read_write(l_lines,name,idp_filter_length)
        prot_len=iup_cand_write[1]
        #gba_cand_write=lcr_extract(lcr_lines,name,gba_filter_length,prot_len)
        
        ##Makes a master span list that itertools.izip_longest uses.
        master_span_list=[iup_cand_write[0],vsl2_cand_write[0],espr_n_cand_write[0],\
                          espr_x_cand_write[0]]#,gba_cand_write[0]]

        ##Prints all info into summary files from all the various dictionaries.
        o = open("summary."+str(protein[0]),"w+")
        print("Protein name:",name,file=o)
        print("Sequence length:",iup_cand_write[1],file=o)
        print("Plaac-score:",plaac_dict[name][0],file=o)
        print("Vitarbi parse:",plaac_dict[name][2],file=o)
        print("Vitarbi centre residue:",plaac_dict[name][1],file=o)
        print("Papa-score:",papa_dict[name][0],file=o)
        print("Papa-score position:",papa_dict[name][1],file=o)
        print("====================================================",file=o)
        print("Filtering value ID:",idp_filter_length,file=o)
        print("Filtering value LC:",gba_filter_length,"\n",file=o)
# Tog bort GBA str
        print("IUPred\tVSL2b\tEspritz N\tEspritz X",file=o)
        print("____________________________________________________",file=o)
        for x in itertools.izip_longest(*master_span_list,fillvalue=""):
                print('\t'.join([str(e) for e in x]),file=o)
        print("====================================================",file=o)
        print("IUPred binary:",iup_cand_write[2],file=o)
	print("VSL2b binary:", vsl2_cand_write[2],file=o)
        print("Espritz N binary:", espr_n_cand_write[2],file=o)
        print("Espritz X binary:", espr_x_cand_write[2],file=o)
        #print("GBA binary:",gba_cand_write[1],file=o)
        o.close()
##Makes a new folder to contain the first candidates summary files.
os.makedirs("first_sweep_summaries_def")
#os.makedirs("second_sweep_summaries_def")


#with open("candidate_prions_first_sweep_def") as first_sweep:
#        #second_sweep_file=open("candidate_prions_second_sweep_def","w+")
#        #print("Systematic Name:\tPlaac LLR-score:\tPapa score",file=second_sweep_file)
#        first_sweep_rows=first_sweep.readlines()
#        temp_sort_list=[]
#        for row in first_sweep_rows[1:]:
#            row=row.split('\t')
#            prot_name=row[0].strip()
#            syst_name=prot_name.split(' ')[0]
#            shutil.move(os.getcwd()+"/summary."+syst_name,\
#                        os.getcwd()+"/first_sweep_summaries_def")
#            if plaac_dict[prot_name][1] in gba_dict[prot_name]:
#                temp_sort_list.append([prot_name,plaac_dict[prot_name][0],papa_dict[prot_name][0]])
#                shutil.move(os.getcwd() + "/first_sweep_summaries_def/summary." + syst_name, \
#                            os.getcwd() + "/second_sweep_summaries_def")
#        temp_sort_list.sort(key=itemgetter(1),reverse=True)
#        for id in temp_sort_list:
#            print(id[0]+"\t"+str(id[1])+"\t"+str(id[2]),file=second_sweep_file)
