from Bio import SeqIO
import re
import sys

###
# Set up inputs / outputs
###
# Read in genbank file for feature information
gb_file = "sequence.gb"
fasta_file = "sequence.fasta" # Note, fasta file should be unannotated
key_file = "key.csv" #key_file = None # Can be set to none to avoid key mapping locus tags to gene names
filter_file = "down.txt"
filtered_out_file = "matching_genes_down.txt" 
motif_match_file = "motif_present.txt"


###
# Read in genome
###

# Read in genome as a single string
fasta_str = ""

with open(fasta_file, "r") as o_fa_file:
    next(o_fa_file)
    for line in o_fa_file:
        fasta_str += line.rstrip()

###
# Parse genbank file for genomic motifs present
# in between CDS features
###

# Genomic motifs to search for
re_list = [
    re.compile(s)
    for s in [
        "TT[AT]A[TAG]A[AGT][TA]TA[GAT]TTAA[TA][TACG]",
        "[TAGC][TA]TTAA[ATCG]TA[ATG][TCGA]T[TAC]T[AT]AA",
    ]
]

# Store motifs with matching genetic feature in inter CDS region
matched_features = set()
for record in SeqIO.parse(gb_file, "genbank"):
    # Currently only expecting one record, but could be exteneded easily
    print("###### RECORD #######")
    print(record)
    print("####################")

    # Essentially implementing a sliding window algorihtm
    prior_stop = None
    current_start = 0
    current_stop = 0

    for feature in record.features:
        if feature.type == "source":  # Skip source features - just overview of geneome
            continue

        if feature.type == "CDS": # Only interface with CDS features
            prior_stop = current_stop
            current_start, current_stop = feature.location.start, feature.location.end

            print("====Feature====")
            print(feature)
            print("CDS fasta sequence")
            print(fasta_str[current_start:current_stop])

            print("prior_start: ", prior_stop)
            print("Current_start: ", current_start)
            print("Current_stop", current_stop)

            # Scan inter CDS region for genomic motifs of interest
            print("inter CDS fasta sequence")
            print(fasta_str[prior_stop:current_start])

            for rgx in re_list:
                matches = rgx.findall(fasta_str[prior_stop:current_start])
                if len(matches) != 0:
                    print("Regex: ", rgx)
                    
                    # Record either gene or locus_tag depending on what's available
                    id_str = "NA"
                    q_dict = feature.qualifiers
                    if "gene" in q_dict:
                        id_str = q_dict["gene"][0]
                    elif "locus_tag" in q_dict:
                        id_str = q_dict["locus_tag"][0]

                    print("MATCH FOUND: ", id_str)
                    matched_features.add(id_str)

# Rename matched motifs according to a key file.
# Should be in format "locus_tag, gname" if trying to replace gene
# name.
if key_file is not None:
    key_dict = {}
    
    with open(key_file) as o_k:
        for line in o_k:
            spl_line = line.rstrip().split(",")

            key_dict[spl_line[0]] = spl_line[1]


    print("Key Dictionary--")
    print(key_dict)
    temp_features = set()

    for f in matched_features:
        if f in key_dict:
            print(f"{f} is actually {key_dict[f]}")
            temp_features.add(key_dict[f])
        else:
            temp_features.add(f)


    matched_features = temp_features


print("========")
print("Features containing genomic motif in inter CDS region")
for f in matched_features:
    print(f)


# Perform feature filtering - only keep features of interest
valid_features = set()
with open(filter_file, "r") as o_file:
    for line in o_file:
        valid_features.add(line.rstrip())


filtered_features = matched_features & valid_features

print("========")
print("Filtered features")
print(filtered_features)

# Record filtered features
with open(filtered_out_file, "w") as o_file:
    if len(filtered_features) == 0:
        o_file.write("No Matches!\n")
    else:
        for f in filtered_features:
            o_file.write(f"{f}\n")


# Record raw features matching motif
with open(motif_match_file, "w") as o_file:
    if len(matched_features) == 0:
        o_file.write("No Matches!\n")
    else:
        for f in matched_features:
            o_file.write(f"{f}\n")
