"""
Created on Jan 22, 2017
@auther: Milinda Samaraweera
This program cleans out the temp files created by pubchem grabber by moving the temp files to the temp_dir
"""
import os
import shutil

def grabber_clean(save_file):

	temp_dir=save_file+"_temp_files"
	
	if os.path.exists(temp_dir):
		shutil.rmtree(temp_dir)
		
	os.makedirs(temp_dir)
	
	try:
		shutil.move("out.txt", os.path.join(temp_dir, 'out.txt') )
	except:
		pass

	try:
		shutil.move(save_file+"_filtered.sdf",os.path.join(temp_dir, save_file+"_filtered.sdf"))
	except:
		pass
		
	try:
		shutil.move(save_file+".sdf.gz",os.path.join(temp_dir, save_file+".sdf.gz"))
	except:
		pass
		
	try:
		shutil.move("molconnIn.sdf",os.path.join(temp_dir,"molconnIn.sdf"))
	except:
		pass
	
	try:
		shutil.move("molconnOut.txt",os.path.join(temp_dir, "molconnOut.txt"))
	except:
		pass
	
	try:
		shutil.move(save_file+"_RI_ECOM50_prefilteredonly.sdf",os.path.join(temp_dir, save_file+"_RI_ECOM50_prefilteredonly.sdf"))
	except:
		pass
		
	try:
		shutil.move(save_file+"_RI_failed.sdf",os.path.join(temp_dir, save_file+"_RI_failed.sdf"))
	except:
		pass
		
	try:
		shutil.move(save_file+"_log.txt",os.path.join(temp_dir, save_file+"_log.txt"))
	except:
		pass
		
	try:
		shutil.move(save_file+"_eliminated_by_RI.sdf",os.path.join(temp_dir, save_file+"_eliminated_by_RI.sdf"))
	except:
		pass
		


def main():
	pass

if __name__ == "__main__":
    main()
			