"""
Created on Jan 22, 2017
@auther: Milinda Samaraweera
This program contains the function to apply the pre-filter on a selected SD file
"""
import shutil
import re
import os
import time
import datetime
import warnings
import gzip
import time
import re
import sys
import requests
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import AllChem

def filter_prefilter(sdf_gz,save_file,high_mass,low_mass,type,log_file):
	
	bio_atoms = ['C', 'H', 'N', 'O', 'P', 'S']
	heavy_isotopes = Chem.MolFromSmarts("[8C,9C,10C,11C,13C,14C,15C,16C,17C,18C,19C,20C,21C,22C,2H,3H,4H,5H,6H,7H,10N,11N,12N,13N,15N,16N,17N,18N,19N,20N,21N,22N,23N,24N,25N,12O,13O,14O,15O,17O,18O,19O,21O,22O,23O,24O,24P,25P,26P,27P,28P,29P,30P,32P,33P,34P,35P,36P,37P,38P,39P,40P,41P,42P,43P,44P,45P,46P,26S,27S,28S,29S,30S,31S,33S,34S,35S,36S,37S,38S,39S,40S,41S,42S,43S,44S,45S,46S,47S,48S,49S]")
	positive_mode= Chem.MolFromSmiles("N")
	alpha_beta_ketone = Chem.MolFromSmiles("C=CC=O")
	ms=[] #Array to store mol objects
		
	if type=='zipped':
		inf = gzip.open(sdf_gz)
		gzsuppl = Chem.ForwardSDMolSupplier(inf)
		
	elif type=='SD':
		gzsuppl  = Chem.SDMolSupplier(sdf_gz)
			
	#print("Candidates removed after mass filter \n Mass\tCandidateID\n")
	log_file.write("\nCandidates removed after mass filter \n Mass\tCandidateID\n")	
	
	for idx,mol in enumerate(gzsuppl):
		if mol is None:
			print ("Error Reading SD entry: " + str(idx))
			log_file.write("Error Reading SD entry " + str(idx)+'\n')
		
		elif mol is not None:
			if round(Descriptors.ExactMolWt(mol),4)>float(low_mass) and round(Descriptors.ExactMolWt(mol),4)<float(high_mass):
				ms.append(mol)
				
			else:
				try:
					pass
					#print(str(round(Descriptors.ExactMolWt(mol),4))+"\t"+mol.GetProp("PUBCHEM_COMPOUND_CID"))
					#log_file.write(str(round(Descriptors.ExactMolWt(mol),4))+"\t"+mol.GetProp("PUBCHEM_COMPOUND_CID")+'\n')
			
				except:
					pass
					#print(str(round(Descriptors.ExactMolWt(mol),4))+"\t"+mol.GetProp("DATABASE_ID"))
					#log_file.write(str(round(Descriptors.ExactMolWt(mol),4))+"\t"+mol.GetProp("DATABASE_ID")+'\n')
					
	if len(ms)==0:
		print("Zero candidates remaning after mass check, exiting...")
		log_file.write("Zero candidates remaning after mass check, exiting...")
		sys.exit()
			
	else:
		print ("number of entries in the updated sdf file after mass check "+save_file+' is: '+str(len(ms)))
		log_file.write("number of entries in the updated sdf file after mass check "+save_file+' is: '+str(len(ms))+'\n')
		
	input_try=0
		
	while input_try<3:
		
		try:
			remove_disconnected = input("Remove candidates having disconnnected structures (e.g. NH4.Cl)? y/n: ")
			remove_charge = input("Remove candidates having an overall charge (e.g. +1,-1,etc.)? y/n: ")
			remove_heavyatoms = input("Remove candidates having heavy isotopes (e.g. 13C, 2H)? y/n: ")
			contain_CHNOPS = input("Retain candidates composed only of CHNOPS? y/n: ")
			keep_pos_ion = input("Retain candidates which can be ionized by ESI in the positive ion mode? y/n: ")
			remove_stereoisomers = input("Remove diastereomers of candidates? y/n: ")
			break
			
		except:
			print("Invalid entry(s)...")
			input_try=input_try+1
				
	if input_try>3:
		print("Inavalid entries, exiting the program...")
	else:
		print("Applying selected prefilters...")
		
	ms_new = []
		
	if remove_disconnected == 'y':
		for mol in ms:
			disc_found = 0
			if str(Chem.MolToSmiles(mol)).find('.') > -1:
				disc_found=1
	
			if disc_found==0:
				ms_new.append(mol)
					
		ms = ms_new
			
		print("Number of Candidates remaning after removel of disconnected structures: "+str(len(ms)))
		log_file.write("\nNumber of Candidates remaning after removel of disconnected structures: "+str(len(ms))+'\n')
		
	ms_new = []
		
	if remove_charge == 'y':
		for mol in ms:
			charge_found = 0
			if rdmolops.GetFormalCharge(mol) != 0:
				charge_found = 1
					
			if charge_found == 0:
				ms_new.append(mol)
					
		ms = ms_new
		print("Number of Candidates remaning after removel of molecules with an overall charge: "+str(len(ms)))
		log_file.write("Number of Candidates remaning after removel of molecules with an overall charge: "+str(len(ms))+'\n')
		
	ms_new = []
		
	if remove_heavyatoms == 'y':
		for mol in ms:
			heavy_isotope_found = 0
				
			if mol.HasSubstructMatch(heavy_isotopes)==True:
				heavy_found=1
				break
					
			try:	
				if int(mol.GetProp("PUBCHEM_ISOTOPIC_ATOM_COUNT"))==0 and heavy_isotope_found==0:
					ms_new.append(mol)
						
			except:
				if heavy_isotope_found==0:
					ms_new.append(mol)
					
		ms = ms_new

		print("Number of Candidates remaning after removel of molecules with heavy isotopes... "+str(len(ms)))
		log_file.write("Number of Candidates remaning after removel of molecules with heavy isotopes... "+str(len(ms))+'\n')
		
	ms_new=[]
		
	if contain_CHNOPS == 'y':
		for mol in ms:
			nonbio_found = 0
			mol_smi = Chem.MolFromSmiles(Chem.MolToSmiles(mol)) 
			for atom in mol_smi.GetAtoms():
				if str(atom.GetSymbol()) not in bio_atoms:
					nonbio_found=1
					break
					
			if nonbio_found==0:
				ms_new.append(mol)
					
		ms = ms_new
		print("Number of Candidates remaning after keeping only with CHNOPS: "+str(len(ms)))
		log_file.write("Number of Candidates remaning after keeping only with CHNOPS: "+str(len(ms))+'\n')
			
	ms_new=[]
		
	if keep_pos_ion  == 'y':
	
		for mol in ms:
			positive_mode_check = 1
		
			if mol.HasSubstructMatch(positive_mode)==False and mol.HasSubstructMatch(alpha_beta_ketone)==False and rdmolops.GetFormalCharge(mol) <=0:
				positive_mode_check=0
					
			if positive_mode_check==1:
				ms_new.append(mol)
					
		ms = ms_new
		print("Number of remaning after keeping candidates detected in the positive ion mode: "+str(len(ms)))
		log_file.write("Number remaning after keeping candidates detected in the positive ion mode: "+str(len(ms))+"\n")
			
	ms_new=[]
							
	if remove_stereoisomers == 'y':
		
		#isomer_out=open(save_file+"_stereoisomer_PCIDs.txt", 'w')
			
		ms_final = []
		smile_list = []
		#isomer_id_list = []
			
			
		for mol in ms:
			#isomer_found = 0
				
			mol_smiles = Chem.MolToSmiles(mol, False, False, -1, True)
				
			if mol_smiles not in smile_list:
				ms_final.append(mol)
				smile_list.append(mol_smiles)
					
			#isomer_id_list.append(mol.GetProp("PUBCHEM_COMPOUND_CID"))	
			#isomer_out.write(str(mol.GetProp("PUBCHEM_COMPOUND_CID")))
				
			#for mol_2 in ms:
					
				#mol_smiles_2 = Chem.MolToSmiles(mol_2, False, False, -1, True)
							
				#if mol_2.GetProp("PUBCHEM_COMPOUND_CID") not in isomer_id_list and mol_smiles==mol_smiles_2:
					#isomer_out.write("\t"+str(mol_2.GetProp("PUBCHEM_COMPOUND_CID"))+"\t")
					#isomer_id_list.append(mol_2.GetProp("PUBCHEM_COMPOUND_CID"))
					#isomer_found=1
								
					#isomer_out.write("\n")
			
		print("Number of Candidates remaning after removing stereoisomers: "+str(len(ms_final)))
		log_file.write("Number of Candidates remaning after removing stereoisomers: "+str(len(ms_final))+'\n')
		#isomer_out.close()
		
	try:
		os.remove(save_file+'_filtered.sdf')
	except: 
		pass

	if remove_stereoisomers=='y':
		print("Resulting candidates saved to "+save_file+'_filtered.sdf '+'and the number of canidates saved was: '+str(len(ms_final)))
			
	else:
		print("Resulting candidates saved to "+save_file+'_filtered.sdf '+'and the number of canidates saved was: '+str(len(ms)))
		
	if remove_stereoisomers=='y':
		w = Chem.SDWriter(save_file+'_filtered.sdf')
		for m in ms_final: w.write(m)
			
	else:
		w = Chem.SDWriter(save_file+'_filtered.sdf')
		for m in ms: w.write(m)
			
	w.close()	
	
def main():
	pass
	
if __name__ == "__main__":
    main()

