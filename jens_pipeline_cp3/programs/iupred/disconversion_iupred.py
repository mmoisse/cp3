import sys
##IUPRED VERSION

##Reads directly from stdin send from the bash script.
##Makes sure PATH error message is not read. Reads all scores as floats.
floats=[]
if __name__ == "__main__":
    for line in sys.stdin:
	if line[0:3]=="IUP":
            continue
        line=line.split("    ")
        line=line[1].strip()
        floats.append(float(line))


    seq_length=0
    binary_map=""
    disorder_residues=0

    ##Read each residue score and create binary map
    for decimal in floats:
        seq_length+=1
        if decimal >=0.5:
            disorder_residues+=1
            binary_map+="D"
        else:
            binary_map+="-"

    spans=[]
    logic_gate_b= "on"
    char_count=0


    ##Read binary map and append start and stop of disorder regions to "spans".
    for char in binary_map:
        char_count+=1
        if char == "D" and logic_gate_b== "on":
            spans.append(char_count)
            logic_gate_b= "off"
        if char == "-" and logic_gate_b== "off":
            spans.append(char_count-1)
            logic_gate_b= "on"
            
    if binary_map[-1]== "D":
        spans.append(char_count)


    ##Read "spans" and for every start of region, print start-stop + length
    residue_counter=0
    for span in spans:
        if residue_counter % 2 == 0:
            print "#" + str(spans[residue_counter]) + "-" + \
                  str(spans[residue_counter+1]) + ":" + \
                  str(spans[residue_counter+1] - spans[residue_counter] + 1)
        residue_counter+=1

    ##Wraps binary map with width 80, prints binary map and sequence length.    
    wrap_count=0
    wrap_list=[]
    wrap_temp=""
    print "====================="
    for x in binary_map:
        wrap_temp+=x
        wrap_count+=1
        if wrap_count % 80 == 0:
            wrap_list.append(wrap_temp)
            wrap_temp=""
    wrap_list.append(wrap_temp)
    
    for y in wrap_list:
        print y
    print "====================="
    print "Sequence length:",seq_length












    
