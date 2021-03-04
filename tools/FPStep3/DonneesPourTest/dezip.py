#!/usr/bin/env python3
# -*-coding:Utf-8 -*
import os
import sys
import gzip

def dezip(file1, out_file1):
 
    #file_out2 = open(out_file2, "w")

    #file_zip = gzip.GzipFile("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv.gz", 'rb')
	file_zip = gzip.GzipFile(file1, 'rb')
	s = file_zip.read()
	file_zip.close()

	#file_dezip = open ("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv", 'wb')
	file_dezip = open (out_file1, 'wb')
	file_dezip.write(s)
	file_dezip.close()

dezip(sys.argv[1], sys.argv[2])

