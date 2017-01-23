#!/usr/bin/env python

"""
This program can be used to download a set of candidates based on exactmass (MIMW) of an unknown within the mass error  
"""

import xml.etree.ElementTree as ET
import xml.etree.ElementTree
from urllib.request import urlopen
import urllib3
import shutil
import re
import os
import time
import datetime
import warnings
import gzip
import time
import sys
import requests
import subprocess
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from prefilter import filter_prefilter
from molconn_v3 import run_molconn, getmolconn_output, ri_filter
from magmainput import magma_plus_input
from cfminput import cfm_input
from pubchemgrabber_cleanup import grabber_clean

warnings.filterwarnings("ignore")

def pug_request(data):
	html_string  = '<?xml version="1.0"?>\n'
	html_string += '<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "https://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">\n'
	html_string += '<PCT-Data>\n'
	html_string += '  <PCT-Data_input>\n' 
	html_string += '    <PCT-InputData>\n'
	html_string += '      <PCT-InputData_download>\n'
	html_string += '        <PCT-Download>\n' 
	html_string += '          <PCT-Download_uids>\n'
	html_string += '            <PCT-QueryUids>\n'
	html_string += '              <PCT-QueryUids_ids>\n'
	html_string += '                <PCT-ID-List>\n'
	html_string += '                  <PCT-ID-List_db>pccompound</PCT-ID-List_db>\n'
	html_string += '                  <PCT-ID-List_uids>\n'  
	for cid in data:
			html_string +='                <PCT-ID-List_uids_E>%s</PCT-ID-List_uids_E>\n' % cid
	html_string += '                  </PCT-ID-List_uids>\n'
	html_string += '                </PCT-ID-List>\n'
	html_string += '              </PCT-QueryUids_ids>\n'
	html_string += '            </PCT-QueryUids>\n'
	html_string += '          </PCT-Download_uids>\n'
	html_string += '          <PCT-Download_format value="sdf"/>\n'
	html_string += '          <PCT-Download_compression value="gzip"/>\n'
	html_string += '        </PCT-Download>\n'
	html_string += '      </PCT-InputData_download>\n'
	html_string += '    </PCT-InputData>\n'
	html_string += '  </PCT-Data_input>\n'
	html_string += '</PCT-Data>\n'
	
	return html_string
	
def pug_query(cmpnd_list):

	xml_obj = xml.etree.ElementTree.fromstring(pug_request(cmpnd_list))
	xml_req = xml.etree.ElementTree.tostring(xml_obj, encoding='us-ascii', method='xml')
	
	r = urlopen('https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi', xml_req)
	out = r.read().decode("utf-8")
	r.close()
	
	while '<PCT-Status value="queued"/>' in out or '<PCT-Status value="running"/>' in out:
		time.sleep(0.5)
		m = re.search("<PCT-Waiting_reqid>(\d+)</PCT-Waiting_reqid>", out)
		if m:
			reqid = m.group(1)
			out = check_query(reqid)
	return out

def check_query(reqid):
	
	check_temp = '''<PCT-Data>
		<PCT-Data_input>
			<PCT-InputData>
					<PCT-InputData_request>
						<PCT-Request>
							<PCT-Request_reqid>%(reqid)s</PCT-Request_reqid>
							<PCT-Request_type value="status"/>
						</PCT-Request>
					</PCT-InputData_request>
				</PCT-InputData>
			</PCT-Data_input>
		</PCT-Data>'''
		
	check = check_temp % {'reqid':str(reqid)}
	
	xml_obj_2 = xml.etree.ElementTree.fromstring(check)
	xml_req_2 = xml.etree.ElementTree.tostring(xml_obj_2, encoding='us-ascii', method='xml')
	
	f = urlopen('https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi',xml_req_2, 2.5)
	out = f.read().decode("utf-8")
	f.close()
	return out

def main():
	
	save_file = str(input("Please enter an identifier to save the work: "))

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
		print("Please enter appropiate values for mass and mass error...")
		log_file.write("Please enter appropiate values for mass and mass error...\n")
		sys.exit()

	low_mass = round(input_mass-mass_diff,4)
	high_mass = round(input_mass+mass_diff,4)

	compound_list = []

	print("Connecting to the PubChem Server to retreive the candidates...")

	try:
		http = urllib3.PoolManager()
		url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&retmax=100000&term='+str(low_mass)+':'+str(high_mass)+'[exactmass]'

		with http.request('GET', url, preload_content=False) as r, open('out.txt', 'wb') as out_file:       
			shutil.copyfileobj(r, out_file)

	except:
		print("Error connecting to PubChem server, please try again later...")
		log_file.write("\nError connecting to PubChem server, please try again later...\n")
		sys.exit()

	r.release_conn()

	with open ('out.txt', 'r') as in_file:
		for line in in_file:
			line = line.lstrip()
			line = line.rstrip()
			if line.startswith('<Id>'):
				line = re.sub('\<Id>', '', line)
				line = re.sub('\</Id>', '', line)
				compound_list.append(int(line))
	in_file.close()
			
	if len(compound_list)==0:
		print("Search resulted in zero canidates, exiting the program...")
		sys.exit()

	
	print("Number of candidates retreived from PubChem matching "+str(input_mass)+" in the masswindow of "+str(input_mass_error_ppm)+" ppm was "+str(len(compound_list)))
	log_file.write("\nNumber of candidates retreived from PubChem matching "+str(input_mass)+" in the masswindow of "+str(input_mass_error_ppm)+" ppm was "+str(len(compound_list))+'\n')

	if len(compound_list)>0:
		result = pug_query(compound_list)
		url = re.search( "<PCT-Download-URL_url>(.*)</PCT-Download-URL_url>", result).group(1)
		print('Downloading sdf files from URL: '+url)
		log_file.write('SD file Downloaded from URL: '+url+'\n')

		try:
			os.remove(save_file+'.sdf.gz')
		except: 
			pass

		sdf = urlopen(url)
		tempfile = open(save_file+".sdf.gz", "wb")
		tempfile.write(sdf.read())
		tempfile.close()
		sdf.close()
		
		filter_prefilter(str(save_file+'.sdf.gz'),save_file,high_mass,low_mass,'zipped',log_file)
		
	else:
		print("Zero candidates resulted from PubChem search, exiting...")
		log_file.write("\nZero candidates resulted from PubChem search, exiting...\n")
		sys.exit()
		
	#This portiion of the code prepares the molconn input file, sends the SD file to metabolomics.pharm.uconn.edu and retreive the predicted RI and Ecom50 values
	
	print("Prepearing MolConn input file...")
	log_file.write("\nPrepearing MolConn input file...\n")

	molconnIn=open('molconnIn.sdf','w')
	suppl2 = Chem.SDMolSupplier(save_file+'_filtered.sdf')

	w_molconn = Chem.SDWriter(molconnIn)
	print(str(len(suppl2))+" Candidates read from the SD saved after adding the pre-filters...")
	log_file.write(str(len(suppl2))+" Candidates read from the SD saved after adding the pre-filters...\n")

	counter = 0

	try:
		os.remove('molconnIn.sdf')
	except:
		pass

	for mol in suppl2:
		if mol is None: 
			continue

		else:
		
			d = mol.GetPropsAsDict()
			mol_smi = str(Chem.MolToSmiles(mol))
			mol_molconn = Chem.MolFromSmiles(mol_smi)
		
			if str(d.get('PUBCHEM_IUPAC_CAS_NAME', None))!='None':
				mol_molconn.SetProp("CompoundName", str(d.get('PUBCHEM_IUPAC_CAS_NAME', None)))
			else:
				#Unable to read IUPAC name from PubChem records using IUPAC InchI instead
				mol_molconn.SetProp("CompoundName", str(d.get('PUBCHEM_IUPAC_INCHI', None)))
				
			mol_molconn.SetProp("Identifier", str(int(d.get('PUBCHEM_COMPOUND_CID', None))))
			mol_molconn.SetProp("InChI", str(d.get('PUBCHEM_IUPAC_INCHI', None)))
			mol_molconn.SetProp("InChIKey", str(d.get('PUBCHEM_IUPAC_INCHIKEY', None)))
			mol_molconn.SetProp("MolecularFormula", str(CalcMolFormula(mol)))
			mol_molconn.SetProp("MonoisotopicMass", str(round(float(Descriptors.ExactMolWt(mol)),4)))
			mol_molconn.SetProp("SMILES", str(Chem.MolToSmiles(mol)))
		
			AllChem.Compute2DCoords(mol_molconn)
		
			w_molconn.write(mol_molconn)
			counter=counter+1

	w_molconn.close()
	molconnIn.close()

	print("Running openbabel to re-format the molconn input file...")
	log_file.write("Running openbabel to re-format the molconn input file...\n")

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

	print("MolConn input file prepared having "+str(after_babel)+" candidates...")
	log_file.write("MolConn input file prepared having "+str(after_babel)+" candidates...\n")

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

	
