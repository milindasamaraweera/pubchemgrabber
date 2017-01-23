"""
Created on Jan 21, 2017
@auther: Milinda Samaraweera
This program creates input files using data from SD files to be used in creating candidate lists to be used in cfm-id
"""
from rdkit import Chem
import csv
import sys

def cfm_input(sd_file,save_file,filter):
	
	f = open(save_file+'_'+filter+'_input_cfmid.txt', 'wt')
	writer = csv.writer(f,delimiter='\t',lineterminator='\n')
	suppl = Chem.SDMolSupplier(sd_file)
	
	if len(suppl)==0:
		print("Zero candidates remaining after RI filter...")
		sys.exit()
	
	for idx,mol in enumerate(suppl):
		
		if mol is None:
				print ("Error Reading SD entry: " + str(idx))
				log_file.write("Error Reading SD entry " + str(idx)+'\n')
				continue
				
		d = mol.GetPropsAsDict()
			
		try:
			Identifier=int(d.get('Identifier', None))
		except:
			Identifier=(str(d.get('Identifier', None)))
			
		writer.writerow((Identifier, str(Chem.MolToSmiles(mol))))
		
	f.close()	

def main():
	pass

if __name__ == "__main__":
    main()
	