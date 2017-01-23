"""
Created on Jan 20, 2017
@auther: Milinda Samaraweera
This program creates csv files using data from SD files to be used in creating .db files for magmaplus
"""
from rdkit import Chem
import csv
import sys

def magma_plus_input(sd_file,save_file,filter):
	
	f = open(save_file+'_'+filter+'_input_magmaplus.csv', 'wt')
	writer = csv.writer(f,quoting=csv.QUOTE_NONNUMERIC,lineterminator='\n')
	suppl = Chem.SDMolSupplier(sd_file)
	writer.writerow( (str("Identifier"), "CompoundName", "MonoisotopicMass","MolecularFormula","SMILES","InChI","InChIKey") )
	
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
			
		writer.writerow((Identifier,str(d.get('CompoundName', None)),round(float(d.get('MonoisotopicMass', None)),4),str(d.get('MolecularFormula', None)),str(d.get('SMILES', None)),str(d.get('InChI', None)),str(d.get('InChIKey', None))))
		
	f.close()	

def main():
	pass

if __name__ == "__main__":
    main()
	