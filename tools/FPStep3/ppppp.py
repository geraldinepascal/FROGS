import sys

file = open(sys.argv[1], "r")

line = file.readline()
while line :
    print("---", line[:5])
    if len(line) and line[0] == "-":
        print(line[:5])
        name_file = line.split(":")[1]
        file_out  = open("v2_" + name_file, "w")

        line = file.readline()
        print(line[:5])
        while line and len(line) and line[0] != "-" :
            print(line[:5])
            file_out.write(line)
            line = file.readline()

        file_out.close()

    if len(line) and line[0] != "-":
        line = file.readline()
    else:
        line = file.readline()
        
file.close()