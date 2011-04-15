from Bio import SeqIO
import sys,os

handle = open(sys.argv[1], "rU")
fw1=open(sys.argv[2],'w')
fw2=open(sys.argv[2]+".len",'w')
for record in SeqIO.parse(handle, "embl") :
    fw1.write(">{}\n{}\n".format(record.id, record.seq))
    fw2.write("{}\t{}\n".format(record.id, len(record.seq)))
handle.close()
fw1.close()
fw2.close()

os.system("cat *.len > all.chr_len")
