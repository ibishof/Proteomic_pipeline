#!/usr/bin/python

import getopt, sys
import string

def badFilename(path):
	return '/' in path

# Set up default parameters
# Optional
mod = "EntrezGene"
permutations = "2000"
method = "z-score"
zscore = "1.96"
pval = "0.05"
num = "3"
# Required
species = None
input = None
denom = None
dataToAnalyze = 'both'
version = 'current' ### Example: EnsMart56Plus

def checkargs():
        ### If a space is present in the denominator or input name, then multiple argmuments for each option will be present
        di = goelite_args.index('--denom'); ii = goelite_args.index('--input')
        if (ii-di) > 2: # Occurs when there is a space in the denominator filename
                goelite_args[di+1] = string.join(goelite_args[di+1:ii]) # Replace the first denominator file string index with a combined index
                del goelite_args[di+2:ii] # delete the leftover parts
        ii = goelite_args.index('--input'); mi = goelite_args.index('--mod')
        if (mi-ii) >2: # Occurs when there is a space in the input filename
                goelite_args[ii+1] = string.join(goelite_args[ii+1:mi]) # Replace the first denominator file string index with a combined index
                del goelite_args[ii+2:mi] # delete the leftover parts                
        
# Command line arguments override default parameters
goelite_args = sys.argv[1:]
checkargs()
opts, args = getopt.getopt(goelite_args, '',
				['species=', 'mod=', 'permutations=', 'method=', 'zscore=',
				 'pval=', 'num=', 'input=', 'denom=','dataToAnalyze=','version='])
for opt, val in opts:
	if opt == "--species":
		species = val
	elif opt == "--mod":
		mod = val
	elif opt == "--permutations":
		permutations = val
	elif opt == "--method":
		method = val
	elif opt == "--zscore":
		zscore = val
	elif opt == "--pval":
		pval = val
	elif opt == "--num":
		num = val
	elif opt == "--input":
		input = val
	elif opt == "--denom":
		denom = val
        elif opt == '--dataToAnalyze':
                dataToAnalyze = val
        elif opt == '--version':
                version = val
if species is None:
	print "Species is missing"
	raise SystemExit, 1
if input is None:
	print "Gene list input file is missing"
	raise SystemExit, 1
if denom is None:
	print "Denominator input file is missing"
	raise SystemExit, 1
if badFilename(input) or badFilename(denom):
	print "Illegal filename"
	raise SystemExit, 1
if len(args) != 0:
	print "Unexpected command line argument"
	print sys.argv[1:]
	raise SystemExit, 1

#print debugging statements
f = open( "/home/socr/c/users2/isaach/goelite.log", "a" )

f.write( "goelite.py input params:" )
f.write( "\n" )
f.write( "[" + species + "]" )
f.write( "[" + mod + "]" )
f.write( "[" + permutations + "]" )
f.write( "[" + method + "]" )
f.write( "[" + zscore + "]" )
f.write( "[" + pval + "]" )
f.write( "[" + num + "]" )
f.write( "[" + input + "]" )
f.write( "[" + denom + "]" )

import os, os.path
def moveToSubdir(filename, dirname):
	sameName = False
	if filename == dirname:
		os.rename(filename, "tmp")
		sameName = True
	os.mkdir(dirname)
	if sameName:
		os.rename("tmp", os.path.join(dirname, filename))
	else:
		os.rename(filename, os.path.join(dirname, filename))
	f.write( "["+dirname + '='+filename+ "]")

try:
	here = os.getcwd()
	there = "/usr/local/projects/cytoscape/GO-Elite_120beta"
	moveToSubdir(input, "input_list")
	moveToSubdir(denom, "denominator")
	args = [ "/usr/bin/python",
			"GO_Elite.py",
			"--species", species,
			"--mod", mod,
			"--permutations", permutations,
			"--method", method,
			"--zscore", zscore,
			"--pval", pval,
			"--num", num,
			"--input", os.path.join(here, "input_list"),
			"--denom", os.path.join(here, "denominator"),
			"--output", here,
                        "--version",version,
                        "--dataToAnalyze",dataToAnalyze,
			]
	os.chdir(there)
	os.execv(args[0], args)
except os.error, e:
	print "Operating system error: %s" % str(e)
	raise SystemExit, 1

f.write( "\n" )
f.close()
