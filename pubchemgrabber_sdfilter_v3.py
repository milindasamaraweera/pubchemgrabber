#!/usr/bin/env python

"""
This program can be used with a SD file downloaded either with PubChem or HMDB to select a set of candidates based on exactmass (MIMW) of an unknown within the mass error
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
import subprocess
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import AllChem
from molconn_v3 import *
from prefilter import filter_prefilter
from magmainput import magma_plus_input
from cfminput import cfm_input
from pubchemgrabber_cleanup import grabber_clean

warnings.filterwarnings("ignore")

def main():

	try:
		input_file = str(sys.argv[1])
		
	except:
		print("Please retry passing the name of the SD file as the first argument!")
		sys.exit()
				
	save_file = str(input("Please input an identifier to save the work: "))
	
	try:
		os.remove(save_file+'_log.txt')
	except:
		pass
			
	log_file=open(save_file+'_log.txt','a')

	try:
		input_mass = round(float(input("Please enter the exact mass of the unknown: ")),4)
		input_mass_error_ppm = round(float(input("Please enter the mass error(ppm): ")),4)
		mass_diff = round((input_mass_error_ppm*input_mass)/(1000000.0),4)
	
	except:
		print("Enter appropiate values...")
		sys.exit()
		
	low_mass = round(input_mass-mass_diff,4)
	high_mass = round(input_mass+mass_diff,4)
	
	print("Candidates will be retreived from SD file matching "+str(input_mass)+" in the masswindow of "+str(input_mass_error_ppm))
	log_file.write("\nCandidates will be retreived from SD file matching "+str(input_mass)+" in the masswindow of "+str(input_mass_error_ppm)+"\n")
		
	stand_babel = str(input("Do you want to run babel to standerdize the SD file? (y/n): "))
	log_file.write("Do you want to run babel to standerdize the SD file? (y/n): \n")
	
	suppl=Chem.SDMolSupplier(input_file)
	stand_before_babel=len(suppl)
	
	if len(suppl)==0:
		print("Zero entries in the SD file, nothing to do, exiting the program...")
		log_file.write("Zero entries in the SD file, nothing to do, exiting the program...\n")
		sys.exit()
		
	else:
		print(str(len(suppl))+" Entries read from the "+input_file)
		log_file.write(str(len(suppl))+" Entries read from the "+input_file+"\n")
		
	if stand_babel=='y' or stand_babel=='Y':
	
		try:
			os.remove('temp.sdf')
			os.remove('babel_'+input_file)
			
		except:
			pass
		
		print('standerdzing the SD file using openbabel...')
		log_file.write("The choice was yes, standerdizing using babel.exe\n")
		
		os.system('babel.exe -isdf '+input_file+' -osdf temp.sdf')
		os.system('rename temp.sdf '+'babel_'+input_file)
		
		suppl=Chem.SDMolSupplier('babel_'+input_file)
		print(str(len(suppl))+" Entries read from the "+input_file)
		log_file.write(str(len(suppl))+" Entries read from the "+'babel_'+input_file+"\n")
		stand_after_babel=len(suppl)
		
		if stand_before_babel==stand_after_babel:
			print("Babel standerdization successfull...")
			log_file.write("Babel standerdization successfull...")
			
		else:
			print("Error in Babel standerdization...")
			log_file.write("Error in Babel standerdization...")
	
	if stand_babel=='y' or stand_babel=='Y':
		filter_prefilter(str('babel_'+input_file),save_file,high_mass,low_mass,'SD',log_file)
		
	else:
		filter_prefilter(input_file,save_file,high_mass,low_mass,'SD',log_file)
	
	#This portiion of the code prepares the molconn input file, sends the SD file to metabolomics.pharm.uconn.edu and retreive the predicted RI and Ecom50 values
	
	print("Prepearing MolConn input file...")
	log_file.write("Prepearing MolConn input file...\n")

	molconnIn=open('molconnIn.sdf','w')
	suppl2 = Chem.SDMolSupplier(save_file+'_filtered.sdf')

	w_molconn = Chem.SDWriter(molconnIn)
	print(str(len(suppl2))+" candidates read frm the SD saved after adding the pre-filters...")
	log_file.write(str(len(suppl2))+" candidates read frm the SD saved after adding the pre-filters...")

	counter = 0

	try:
		os.remove('molconnIn.sdf')
	except:
		pass

	for idx,mol in enumerate(suppl2):
		if mol is None:
			print ("No molecule: " + str(idx))
			continue
		
		else:
	
			d = mol.GetPropsAsDict()
			mol_smi = str(Chem.MolToSmiles(mol))
			mol_molconn = Chem.MolFromSmiles(mol_smi)
			
			if 'JCHEM_IUPAC' in d:
				if str(d.get('JCHEM_IUPAC', None))!='None':
					mol_molconn.SetProp("CompoundName", str(d.get('JCHEM_IUPAC', None)))
				else:
					mol_molconn.SetProp("CompoundName", str(d.get('JCHEM_TRADITIONAL_IUPAC', None)))
				mol_molconn.SetProp("Identifier", str(d.get('HMDB_ID', None)))
				mol_molconn.SetProp("InChI", str(d.get('INCHI_IDENTIFIER', None)))
				mol_molconn.SetProp("InChIKey", str(d.get('INCHI_KEY', None)))
				try:
					mol_molconn.SetProp("MolecularFormula", str(CalcMolFormula(mol)))
				except:
					mol_molconn.SetProp("MolecularFormula", str(d.get('FORMULA', None)))
				try:
					mol_molconn.SetProp("MonoisotopicMass", str(round(float(Descriptors.ExactMolWt(mol)),4)))
				except:
					mol_molconn.SetProp("MonoisotopicMass", str(d.get('EXACT_MASS', None)))
				mol_molconn.SetProp("SMILES", str(Chem.MolToSmiles(mol)))
				
		
			elif 'PUBCHEM_candidate_CID' in d:
		
				if str(d.get('PUBCHEM_IUPAC_CAS_NAME', None))!='None':
					mol_molconn.SetProp("CompoundName", str(d.get('PUBCHEM_IUPAC_CAS_NAME', None)))
				else:
					mol_molconn.SetProp("CompoundName", str(d.get('PUBCHEM_IUPAC_INCHI', None)))
				mol_molconn.SetProp("Identifier", str(int(d.get('PUBCHEM_candidate_CID', None))))
				mol_molconn.SetProp("InChI", str(d.get('PUBCHEM_IUPAC_INCHI', None)))
				mol_molconn.SetProp("InChIKey", str(d.get('PUBCHEM_IUPAC_INCHIKEY', None)))
				try:
					mol_molconn.SetProp("MolecularFormula", str(CalcMolFormula(mol)))
				except:
					mol_molconn.SetProp("MolecularFormula", str(d.get('PUBCHEM_MOLECULAR_FORMULA', None)))
				try:
					mol_molconn.SetProp("MonoisotopicMass", str(round(float(Descriptors.ExactMolWt(mol)),4)))
				except:
					mol_molconn.SetProp("MonoisotopicMass", str(d.get('PUBCHEM_MONOISOTOPIC_WEIGHT', None)))
				mol_molconn.SetProp("SMILES", str(Chem.MolToSmiles(mol)))
				
			elif 'CompoundName' in d:
				if str(d.get('CompoundName', None))!='None':
					mol_molconn.SetProp("CompoundName", str(d.get('CompoundName', None)))
				else:
					mol_molconn.SetProp("CompoundName", str(d.get('CompoundName', None)))
				try:
					mol_molconn.SetProp("Identifier", str(int(d.get('Identifier', None))))
				except:	
					mol_molconn.SetProp("Identifier", str(d.get('Identifier', None)))
				mol_molconn.SetProp("InChI", str(d.get('InChI', None)))
				mol_molconn.SetProp("InChIKey", str(d.get('InChIKey', None)))
				try:
					mol_molconn.SetProp("MolecularFormula", str(CalcMolFormula(mol)))
				except:
					mol_molconn.SetProp("MolecularFormula", str(d.get('MolecularFormula', None)))
				try:
					mol_molconn.SetProp("MonoisotopicMass", str(round(float(Descriptors.ExactMolWt(mol)),4)))
				except:
					mol_molconn.SetProp("MonoisotopicMass", str(d.get('MonoisotopicMass', None)))
				mol_molconn.SetProp("SMILES", str(Chem.MolToSmiles(mol)))
				
		
			AllChem.Compute2DCoords(mol_molconn)
			w_molconn.write(mol_molconn)
			counter=counter+1
		
	w_molconn.close()
	molconnIn.close()

	before_babel=counter

	while True:
		try:
			with open('molconnIn.sdf', 'rb') as _:
				break
		except:
			time.sleep(3)		
	try:
		os.remove('temp.sdf')
	except:
		pass
	
	os.system('babel.exe -isdf molconnIn.sdf -osdf temp.sdf')

	try:
		os.remove('molconnIn.sdf')
	except:
		pass

	os.system('rename temp.sdf molconnIn.sdf')

	while True:
		try:
			with open('molconnIn.sdf', 'rb') as _:
				break
		except: 
			time.sleep(3)	

	after_babel_suppl=Chem.SDMolSupplier('molconnIn.sdf')
	after_babel=len(after_babel_suppl)

	if before_babel==after_babel:
		print("openbabel conversion successfull...")
		log_file.write("openbabel conversion successfull...\n")
	else:
		print("Entries missing after babel conversion...")
		log_file.write("Entries missing after babel conversion...\n")
		sys.exit()

	print("MolConn input file prepared having "+str(after_babel)+" compunds...")
	log_file.write("MolConn input file prepared having "+str(after_babel)+" compunds...")

	if after_babel==0 or counter==0:
		print("Nothing in the molconnIn.sdf terminating the program...")
		log_file.write("Nothing in the molconnIn.sdf terminating the program...\n")
		sys.exit()

	ri_choice=str(input("Do you wish to apply the RI filter? (y/n): "))
	
	if ri_choice=='y' or ri_choice=='Y':
		shutil.copyfile('molconnIn.sdf',save_file+'_prefiltered.sdf')
		run_molconn(save_file,log_file)
		getmolconn_output(save_file,log_file)
		ri_filter(save_file,log_file)
		print("Selected Pre-filters are applied and saved to file: "+save_file+'_prefiltered.sdf')
		log_file.write("Selected Pre-filters are applied and saved to file: "+save_file+'_prefiltered.sdf\n')
		magma_plus_input(save_file+'_prefiltered.sdf',save_file,'PF')
		cfm_input(save_file+'_prefiltered.sdf',save_file,'PF')
		print("Selected Pre-filters and RI filters are applied and saved to file: "+save_file+'_prefiltered_RIfiltered.sdf')
		log_file.write("Selected Pre-filters and RI filters are applied and saved to file: "+save_file+'_prefiltered_RIfiltered.sdf\n')
		magma_plus_input(save_file+'_prefiltered_RIfiltered.sdf',save_file,'PFRI')
		cfm_input(save_file+'_prefiltered_RIfiltered.sdf',save_file,'PFRI')
		
	else:
		shutil.copyfile('molconnIn.sdf',save_file+'_prefiltered.sdf')
		print("Selected Pre-filters are applied and saved to file: "+save_file+'_prefiltered.sdf')
		log_file.write("Selected Pre-filters are applied and saved to file: "+save_file+'_prefiltered.sdf\n')
		magma_plus_input(save_file+'_prefiltered.sdf',save_file,'PF')
		cfm_input(save_file+'_prefiltered.sdf',save_file,'PF')
		
	log_file.flush()
	log_file.close()
	
	grabber_clean(save_file)
	
if __name__ == "__main__":
    main()

