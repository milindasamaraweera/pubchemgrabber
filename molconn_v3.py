#!/usr/bin/env python

import paramiko
import requests
import os
import sys
import time
import string
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def run_molconn(save_file,log_file):
	#Login information to the server as molfind
	username = 'molfind'
	password='molconn@uconn'
	hostname ='metabolomics.pharm.uconn.edu'
	port=22
	
	try:
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		
		print("Connecting to: metabolomics.pharm.uconn.edu")
		log_file.write("\nConnecting to: metabolomics.pharm.uconn.edu")
		
		ssh.connect(hostname, username=username, password=password)
		ssh.get_transport().window_size = 3 * 1024 * 1024
		time.sleep(10)
		
		print("Successfully connected.")
		log_file.write("Successfully connected.\n")
	
	except:
		print("Connection failed, Server might be busy, please try again later...")
		log_file.write("Connection failed, Server might be busy, please try again later...\n")
		
	suppl = Chem.SDMolSupplier('molconnIn.sdf')
	print(str(len(suppl))+" Molecules read from the molconnIn.sdf file...\n")
	log_file.write(str(len(suppl))+" Molecules read from the molconnIn.sdf file...\n")
	
	if len(suppl)==0:
		print("Terminating, zero compounds found in the SD file...")
		log_file.write("Terminating, zero compounds found in the SD file...\n")
		sys.exit()
	
	print("Calculating RI+Ecom50 values on the remote server...")
	log_file.write("Calculating RI+Ecom50 values on the remote server...\n")
	
	sftp = ssh.open_sftp()
	
	try:
		#This will remove any previously created molconn input and output files
		os.remove('molconnOut.txt')
		stdin, stdout, stderr = ssh.exec_command("cd /home/molfind/00_MolFind_predict/pubchemGrabber; rm molconnIn.sdf")
		stdin, stdout, stderr = ssh.exec_command("cd /home/molfind/00_MolFind_predict/pubchemGrabber; rm molconnOut.txt")
			
	except:
		pass
			
	transfer_file_try=0
	
	while transfer_file_try < 3:
		try:
			sftp.put('molconnIn.sdf','/home/molfind/00_MolFind_predict/pubchemGrabber/molconnIn.sdf')
			#stdin, stdout, stderr = ssh.exec_command("cd /home/molfind/00_MolFind_predict/pubchemGrabber/; ls")
			break
			
		except:
			transfer_file_try=transfer_file_try+1
			
		if transfer_file_try>3:
			sys.exit()
			
	if transfer_file_try>3:
		print("Error transfering the molconnIn.sdf file...")
		log_file.write("Error transfering the molconnIn.sdf file...\n")
		sys.exit()
			
	try:
		stdin, stdout, stderr = ssh.exec_command("cd /home/molfind/00_MolFind_predict/; ./04_MolFind_predict.sh /home/molfind/00_MolFind_predict/pubchemGrabber/ Identifier verbose	")
		print("********************************running molconn********************************")
		log_file.write("********************************running molconn********************************\n")
		for line in stdout.read().splitlines():
			print(line)
			log_file.write(str(line)+'\n')
				
		print("********************************molconnfinished********************************")
		log_file.write("********************************molconnfinished********************************\n")
			
	except:			
		print("Error occured while running molconn...")
		log_file.write("Error occured while running molconn...\n")
		sys.exit()
			
	ssh.close()
	sftp.close()
	
def getmolconn_output(save_file,log_file):
	trys=0
	try:
		username = 'molfind'
		password='molconn@uconn'
		hostname ='metabolomics.pharm.uconn.edu'
		port=22
		localpath =os.getcwd()
		t = paramiko.Transport((hostname, port))

		t.connect(username=username, password=password)
		sftp = paramiko.SFTPClient.from_transport(t)
		
	except:
		print("Server busy please try again later...")
		log_file.write("Server busy please try again later...\n")
		
	dirlist = sftp.listdir('/home/molfind/00_MolFind_predict/pubchemGrabber/')

	while 'molconnOut.txt' not in dirlist:
		time.sleep(1)
		trys=trys+1
		dirlist = sftp.listdir('/home/molfind/00_MolFind_predict/pubchemGrabber/')
		
		if trys>900:
			print("Error calculating the RI values, process taking too much time...")
			log_file.write("Error calculating the RI values, process taking too much time...\n")
			sys.exit()
	
	if 'molconnOut.txt' in dirlist:
		print("Transfering predicted Ecom50 and RI values from the server..")
		log_file.write("Transfering predicted Ecom50 and RI values from the server...\n")
		sftp.get('/home/molfind/00_MolFind_predict/pubchemGrabber/molconnOut.txt', localpath +'\\'+'molconnOut.txt')
		print("Transfer Complete...")
		log_file.write("Transfer Complete...\n")
			
	else:
		print("Error in transfering molconnOut.txt")
		log_file.write("Error in transfering molconnOut.txt\n")
		
	t.close()
		
def ri_filter(save_file,log_file):
	
	valid_RI_input='Null'
	
	while valid_RI_input!='True':
		try:
			experimental_RI = float(input("Please enter the experimental RI value: "))
			valid_RI_input='True'
		except:
			print("Invalid input please try again...")
			valid_RI_input='False'
	
	molconn_out = open('molconnOut.txt', 'r')
	
	molconnres=open(save_file+'_RI_ECOM50_prefilteredonly.sdf','w')
	molconn_final=open(save_file+'_prefiltered_RIfiltered.sdf','w')
	molconn_RI_rejected=open(save_file+'_eliminated_by_RI.sdf','w')
	molconn_RI_failed=open(save_file+'_RI_failed.sdf','w')
	
	
	w_molconn = Chem.SDWriter(molconnres)
	final_molconn =Chem.SDWriter(molconn_final)
	RI_rejected = Chem.SDWriter(molconn_RI_rejected)
	RI_failed = Chem.SDWriter(molconn_RI_failed)
	
	counter_annotated=0
	counter_final=0
	
	RI_error = 101 # Enter the error of the model
	suppl3 = Chem.SDMolSupplier('molconnIn.sdf')
	
	print(str(len(suppl3))+" Molecules read from the molconnIn.sdf file...")
	log_file.write("\n"+str(len(suppl3))+" Molecules read from the molconnIn.sdf file...\n")
	
	print('Filtering compounds RI values greater or less than '+str(RI_error)+' RI units the of the experimental value of '+str(experimental_RI))
	log_file.write('Filtering compounds RI values greater or less than '+str(RI_error)+' RI units the of the experimental value of '+str(experimental_RI)+'\n')
	
	for idx,mol in enumerate(suppl3):
		if mol is None:
			print ("Error Reading SD entry: " + str(idx))
			log_file.write("Error Reading SD entry " + str(idx)+'\n')
			continue
			
		else:
			d = mol.GetPropsAsDict()
			mol_smi = str(Chem.MolToSmiles(mol))
			mol_molconn = Chem.MolFromSmiles(mol_smi)
			
			try:
				chem_name_molconnIn = 'NAME;ID_'+str(int(d.get('Identifier', None)))
				
			except:
				chem_name_molconnIn = 'NAME;ID_'+str(d.get('Identifier', None))
			
			compound_found=0
			RI_calculated=0
			
			for line in molconn_out:
				if line.startswith('NAME;ID_'):
					chem_name_molconnOut = line.split()[0]
					
					if chem_name_molconnIn==chem_name_molconnOut:
						compound_found=1
						
						mol_molconn.SetProp("CompoundName", str(d.get('CompoundName', None)))
						try:
							mol_molconn.SetProp("Identifier", str(int(d.get('Identifier', None))))
						except:
							mol_molconn.SetProp("Identifier", str(d.get('Identifier', None)))
						
						mol_molconn.SetProp("InChI", str(d.get('InChI', None)))		
						mol_molconn.SetProp("Predicted_Retention_Index", str(line.split()[2]))
						mol_molconn.SetProp("Predicted_ECOM50", str(line.split()[1]))
						mol_molconn.SetProp("InChIKey", str(d.get('InChIKey', None)))
						mol_molconn.SetProp("MolecularFormula", str(d.get('MolecularFormula', None)))
						mol_molconn.SetProp("MonoisotopicMass", str(round(float(d.get('MonoisotopicMass', None)),4)))
						mol_molconn.SetProp("SMILES", str(d.get('SMILES', None)))
						
						AllChem.Compute2DCoords(mol_molconn)
						
						w_molconn.write(mol_molconn)
						counter_annotated=counter_annotated+1
						
						if str(line.split()[2])!='failed':
							if float(line.split()[2])>experimental_RI-RI_error and float(line.split()[2])<experimental_RI+RI_error:
								final_molconn.write(mol_molconn)
								counter_final=counter_final+1
								RI_calculated=1
								
						elif str(line.split()[2])=='failed':
							final_molconn.write(mol_molconn)
							RI_failed.write(mol_molconn)
							counter_final=counter_final+1
							RI_calculated=1
							
						if RI_calculated==0:
							RI_rejected.write(mol_molconn)
									
			if compound_found==0:
				mol_molconn.SetProp("CompoundName", str(d.get('CompoundName', None)))
				try:
					mol_molconn.SetProp("Identifier", str(int(d.get('Identifier', None))))
				except:
					mol_molconn.SetProp("Identifier", str(d.get('Identifier', None)))
					
				mol_molconn.SetProp("InChI", str(d.get('InChI', None)))		
				mol_molconn.SetProp("Predicted_Retention_Index", str(line.split()[2]))
				mol_molconn.SetProp("Predicted_ECOM50", str(line.split()[1]))
				mol_molconn.SetProp("InChIKey", str(d.get('InChIKey', None)))
				mol_molconn.SetProp("MolecularFormula", str(d.get('MolecularFormula', None)))
				mol_molconn.SetProp("MonoisotopicMass", str(round(float(d.get('MonoisotopicMass', None)),4)))
				mol_molconn.SetProp("SMILES", str(d.get('SMILES', None)))
				
				AllChem.Compute2DCoords(mol_molconn)
				w_molconn.write(mol_molconn)
				counter_annotated=counter_annotated+1
							
			molconn_out.seek(0)
	
	print("MolConn results annotated for "+str(counter_annotated)+" compunds...")
	log_file.write("MolConn results annotated for "+str(counter_annotated)+" compunds...\n")
	print("RI filter resulted in "+str(counter_final)+" compounds...")
	log_file.write("RI filter resulted in "+str(counter_final)+" compounds...\n")
	
	print("RI filtered compounds are appended to :"+save_file+'_prefiltered_RIfiltered.sdf...')
	log_file.write("RI filtered compounds are appended to :"+save_file+'_prefiltered_RIfiltered.sdf...\n')
	
	w_molconn.close()
	final_molconn.close()
	RI_rejected.close()
	RI_failed.close()
	molconnres.close()
	molconn_final.close()
	
def main():
	pass
	
if __name__ == "__main__":
    main()