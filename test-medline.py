from Bio import Entrez
Entrez.email = "jmend3z@gmail.com"  # Always tell NCBI who you are
handle = Entrez.egquery(term="biopython")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    print(row["DbName"], row["Count"])

pmid = "19304878"
record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))
print(record)
"""from Bio import Entrez
handle = open("einfo3.xml", "rb")
record = Entrez.read(handle, validate=False)
handle.close()

from Bio import Medline
with open("Medline/pubmed_result2.txt") as handle:
    records = Medline.parse(handle)
    for record in records:
        print(record['TI'])
"""
