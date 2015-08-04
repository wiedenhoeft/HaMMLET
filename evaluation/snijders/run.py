import os
hammlet = "python ../../hammlet.py"

for line in file("info.txt"):
	line=line.strip()
	gsm, cellLine, repo = line.split("\t")
	callstring = "%s  -a 6 0.1 0.85 1 1 1 -p -L1 -T -b 900 -i 100 model.txt %s.csv -o %s_%s_%s  -x 'Probe number' -y 'Log2 ratio'" % (hammlet, gsm, repo, cellLine, gsm)
	print callstring
	os.system(callstring)
