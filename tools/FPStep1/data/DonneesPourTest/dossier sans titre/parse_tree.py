import re

#trees = loads('(A,B,(C,D)E)F;')
file  = open("out.tree", "r")
#line  = file.readline()
#trees = loads(line)

line = file.readline()
print(line.split(","))
list_cluster = []
while line:
	for i, v in enumerate(line.split(",")):
		group = re.search("(Cluster_[0-9]+)", v)
		if group:
			ide = group.group(1)
			list_cluster.append(ide)
	line = file.readline()
print(list_cluster)
file.close()

file_fasta = open("test2.fasta", "r")
file_out  = open("out.fasta", "w")

line = file_fasta.readline()
while line:
	if line[0] == ">":
		print(line[1:].strip())
		if line[1:].strip() not in list_cluster:
			print(line[1:].strip())
			file_out.write(line[1:].strip()+"\n")

	line = file_fasta.readline()

file_fasta.close()
file_out.close()
