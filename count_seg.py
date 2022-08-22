#Counting the segments for each species and making summary, splitting the fasta files
import sys #Library for getting argument from command line
#Biopython for importing and manipulating sequence records
from Bio import SeqIO
from Bio.Seq import Seq
#Accept and assign the arguments from command line
family=sys.argv[1] #Family name that is being counted
filename=sys.argv[2] #Name of the file with all the genbank fasta sequences
genbank_filename=f"../results/{family}/{family}.csv" #The csv file that contains the metadata for the family
out_file=f"../results/{family}/{family}.txt" #The file with counts for each species and family

species = [] #stores all species name
segments = [] #stores all the segment names`

#Open csv and go through the different lines, each line represents a sequence
#Create lists that contain all of the species and segment names
with open(genbank_filename) as file_handler:
    for line in file_handler:
        line = line.strip() #Removes trailing and leading spaces
        splits = line.split(",") #CSV has to be split by comma, creates a List containing the values of each column
        species_name = splits[1] #Species name in the second column
        #Add species to list
        if species_name not in species:
            species.append(species_name) #Add species name to the list
        segment_name = (((splits[6].replace(" ","")).replace("-","")).replace("_","")).lower() #To avoid redundancy eg. RNA_1, RNA 1, RNA1 are all the same.
        #Some segments have been inconcistenly named
        if segment_name=="small":
            segment_name="s"
        if segment_name=="large":
            segment_name="l"
        if segment_name=="medium" or segment_name=="middle":
            segment_name="m"
        #Add segments to the list if not already existing
        if segment_name not in segments:
            segments.append(segment_name)

data = {} #Dict with the count data of species and segment
isolates = {} #Dictionary with the acc_idss of species and segments
#Creates dictionary within a dictionary
for spe_name in species:
    data[spe_name] = {}
    isolates[spe_name] = {}
    for seg_name in segments:
        data[spe_name][seg_name] = 0 #Count for all initialised as 0 so that there are no empty records, will aid csv
        isolates[spe_name][seg_name] = [] #Empty lists, that will contain acc id for relevant records
    data[spe_name]["total"] = 0 #Will contain total count of sequences for each species
#Counts each segment for each species and saves isolate name
with open(genbank_filename) as file_handler:
    for line in file_handler:
        line = line.strip()
        splits = line.split(",")
        species_name = splits[1]
        segment_name = (((splits[6].replace(" ","")).replace("-","")).replace("_","")).lower()
        if segment_name=="small":
            segment_name="s"
        if segment_name=="large":
            segment_name="l"
        if segment_name=="medium" or segment_name=="middle":
            segment_name="m"
        acc_id = splits[0] #Gets accession id of the sequence
        data[species_name][segment_name] += 1 #Count for particular species AND segment increases by one
        data[species_name]["total"]+=1 #Total count for soecies increase by one
        isolates[species_name][segment_name].append(acc_id) #Adds the accesion id to the list
#Write out the fasta records to separate output files and save filename to a list
with open(out_file,"w") as outwrite:
    for sp in data:
        spe=sp.replace(" ","-").replace("/","-") #Removes spaces as file names cannot have them
        for seg in data[sp]:
            if seg=="total": break
            records=[] #List that will contain the record for given species and segement, that has to be written out to the
            if data[sp][seg]>0:
                for seq_record in SeqIO.parse(filename, "fasta"):
                    #If accession id matches the ones in the list for the species and segment, it will be added to the records
                    if seq_record.id in isolates[sp][seg]:
                        records.append(seq_record)
                #Add all records to a fasta file with the name containing species and segment name seperated by _
                SeqIO.write(records,f"../results/{family}/fastafiles/{spe}_{seg}.fasta","fasta")
                outwrite.write(f"{spe}_{seg}\n") #Save the filename to a list for iterrating
output_filename=f"../results/{family}/{family}_segments.txt" #TSV having counts for the segments of each species in the family before the analysis
#Writes out  the count of each segment to the output file
with open(output_filename, "w") as file_output:
     #Headers for the file
     out_str = "species"
     for seg in segments:
         out_str += "\t" + seg
     file_output.write(out_str + "\tTotal\n")
     #Adding species name and the segment wise counts
     for sp in data:
         out_str = sp
         for seg in data[sp]:
             out_str += "\t" + str(data[sp][seg])
         file_output.write(out_str + "\n")
#Creates the text file in which each species segment count would be entered after filtration step
with open(f"../results/{family}/{family}_filtered.txt","w") as filt:
    filt.write("Species\tSegment\tCount\n")
