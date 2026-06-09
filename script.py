from Bio import SeqIO
import re
import sys


# Read in genbank file for feature information
gb_file = "sequence.gb"
fasta_file = "sequence.fasta"
re_list = [
    re.compile(s)
    for s in [
        "TT[AT]A[TA]A[AGT][TA]TA[GAT]TTAA[TA][TC]",
        "[GA][TA]TTAA[ATC]TA[AT][TCA]T[TA]T[AT]AA",
    ]
]


fasta_str = ""

with open(fasta_file, "r") as o_fa_file:
    next(o_fa_file)
    for line in o_fa_file:
        fasta_str += line.rstrip()

print(fasta_str[0:5])
print(fasta_str[-5:])

print(len(fasta_str))
i = 0

for record in SeqIO.parse(gb_file, "genbank"):
    # Currently only expecting one record, but could likely be exteneded easily
    print("###### RECORD #######")
    print(record)
    print("####################")

    prior_stop = None
    current_start = None
    current_stop = 0

    for feature in record.features:
        if feature.type == "source":  # Skip source features - just overview of geneome
            continue

        if feature.type == "CDS":
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
                print(rgx)
                print(matches)

                if len(matches) != 0:
                    id_str = "NA"
                    q_dict = feature.qualifiers
                    if "gene" in q_dict:
                        id_str = q_dict["gene"]
                    elif "locus_tag" in q_dict:
                        id_str = q_dict["locus_tag"]
                    # print(feature.qualifiers)
                    print("MATCH FOUND ", id_str)

        i += 1
        # if i == 10:
        #    sys.exit()


# with open(gb_file, "r") as o_gb_file:
#    for feature in SeqIO.parse(o_gb_file, "gb").features:
#        print(feature)


# Iterate feature by feature, create an upstream window

# Extract upstream window - compare feature against
