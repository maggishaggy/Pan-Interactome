import os, json, pandas as pd


cwd = os.getcwd()

for file in os.listdir(cwd):
	if file.endswith(".json"):
		json_file = file


info_dict = dict()

with open(json_file) as data:
	cluster_dict = json.load(data)
	
	for cluster in cluster_dict:
		
		for key, value in cluster.iteritems():
			
			if key == "GName" and value.encode('ascii', 'ignore') not in info_dict:
				gene_name = value.encode('ascii', 'ignore')
			if key == "locus":
				strains = value.encode('ascii', 'ignore').split(" ")
		
		info_dict[gene_name] = strains 


strain_set = set()

for value in info_dict.itervalues():

	for item in value:
	
		strain_set.add('_'.join(item.split('_', 2)[:2]))


	
df = pd.DataFrame(0, index=info_dict.keys(), columns=strain_set)

for strain in list(df):

	for gene in df.index:
	
		if strain in " ".join(info_dict[gene]):
			df.at[gene, strain] = 1




