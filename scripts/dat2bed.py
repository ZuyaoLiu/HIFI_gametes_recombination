import json
from argparse import (ArgumentParser, FileType)

def parse_args():
	parser = ArgumentParser(description='Convert the Tandem Repeat Finder (TRF) DAT file into BED/JSON format for genotyping STRs')
	parser.add_argument('--dat', type=str, required=True, help='Input DAT file produced by Tandem Repeat Finder (TRF) using the -d option')
	parser.add_argument('--out', type=str, required=False, help='Output BED/JSON file based on Tandem Repeat Finder (TRF) data')
	parser.add_argument('--tool', type=str, required=True, choices=['expansionhunter', 'lobstr', 'gangstr', 'hipstr', 'repeatseq', 'gatk'], help='Name of the tool for which the bed/json file will be generated')
	
	return parser.parse_args() 

def main():
	args = parse_args()
	datfile = args.dat
	outfile = args.out
	toolname = args.tool

	with open(outfile, 'w') as out:
		chrom = ""
		catalogue = []

		with open(datfile, 'r') as dat:
			for line in dat:
				splitline = line.split()
				if line.startswith("@N"):
					chrom = line.split("@")[1].strip()
				else:
					try:
						try:
							int(splitline[0])
						except ValueError:
							continue
						start = splitline[0]
						end = splitline[1]
						motif_length = splitline[2]
						reference_length = splitline[3]
						motif = splitline[13]
						alignment_score = splitline[7]
						none = '.'

						if toolname == 'expansionhunter':
							catalogue.append({
								"LocusId": chrom + '-' + str(start) + '-' + str(end),
								"LocusStructure": '(' + motif + ')*',
								"ReferenceRegion": chrom + ':' + str(start) + '-' + str(end),
								"VariantType": "Repeat",
								})

						elif toolname == 'lobstr':
							out.write('\t'.join([chrom, start, end, motif_length, reference_length, none, none, none, alignment_score, none, none, none, none, none, motif]) + '\n') # for LobSTR

						elif toolname == 'gangstr':
							out.write('\t'.join([chrom, start, end, motif_length, motif]) + '\n') # for GangSTR
		
						elif toolname == 'hipstr':
							out.write('\t'.join([chrom, start, end, motif_length, reference_length]) + '\n') # for HipSTR

						elif toolname == 'repeatseq':
							out.write(chrom + ':' + start + '-' + end + '\t' + motif_length + '_'.join([reference_length, motif_length, start, end, alignment_score, splitline[8], splitline[9], splitline[10], splitline[11], splitline[12], motif]) + '\n') # for RepeatSeq

						elif toolname == 'gatk':
							out.write(chrom + '\t' + start + '\t' + end + '\n') # for GATK

						else:
							print("Tool not specified")
					except IndexError:
						pass

		if toolname == 'expansionhunter':
			out.write(json.dumps(catalogue, indent = 4)) # for ExpansionHunter


if __name__ == '__main__':
	main()

	print("Conversion finished")

