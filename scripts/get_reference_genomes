from fileinput import filename
from ftplib import FTP
import random
import os
import os.path
import gzip
import shutil
from numpy import save

# login to ftp and direct to needed directory
ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login(user = "anonymous", passwd = "ohr5003") # login creds to the ftp server
ftp.cwd('genomes/.vol2/refseq/bacteria/') 
ftp.encoding = "utf-8"

def goToDirect(): # go to the directory with the FASTA files
    final_Directory = False
    while final_Directory != True:
        for thing in ftp.nlst():
            if thing[-6:] == "fna.gz":
                final_Directory = True
        if final_Directory != True:
            ftp.cwd(ftp.nlst()[0])
    return final_Directory

def downloadFNAfile(path, current_directory_name): # download all FASTA files in a directory
    path = makeFile(path + f"/{current_directory_name}")
    for filename in ftp.nlst():
        if filename[-6:] == "fna.gz":
            completeName = os.path.join(path, filename)         
            file = open(completeName, "wb")
            ftp.retrbinary(f"RETR {filename}", file.write)
            file.close()

def returnToOriginalDirect(): # return to main directory
    ftp.cwd('/') 
    ftp.cwd('genomes/.vol2/refseq/bacteria/') 

def makeFile(save_path): # make folders
    path = save_path 
    try: 
        os.mkdir(path) 
    except OSError as error: 
        print(error) 
    return path

def unZip(path):
    for bacteriaDirectory in os.listdir(path):
        new_path = path + f"/{bacteriaDirectory}"
        for fnaFILE in os.listdir(new_path):
            final_path = new_path + f"/{fnaFILE}"
            with gzip.open(final_path, 'rb') as f_in:
                with open(final_path[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

def removeZip(path):    
    for bacteriaDirectory in os.listdir(path):
        new_path = path + f"/{bacteriaDirectory}"
        for fnaFILE in os.listdir(new_path):
            name, extension = os.path.splitext(fnaFILE)
            if extension == ".gz":
                os.remove(os.path.join(new_path, fnaFILE))

def main(number_of_genemoes, save_path):
    print("Connecting...")
    path = makeFile(save_path + 'refrence_genemoes') # make main folder
    bacteriaDirectoryNames = ftp.nlst() # list of possible refrence bacterias
    randomlist = random.sample(range(0, 41433), number_of_genemoes) # get a random sample
    print("Downloading...")
    for aRandomNumber in randomlist: # main loop to save file
        current_directory_name = bacteriaDirectoryNames[aRandomNumber] # whats the bacteria name
        ftp.cwd(current_directory_name) # go to the specfied bacteria directory
        goToDirect() # go to the directory with the FASTA files
        downloadFNAfile(path, current_directory_name) # download the FASTA files
        returnToOriginalDirect() # return to the main directory once done
    ftp.close() # close the connection
    print("Unzipping files...")
    unZip(path) # unzip the files
    removeZip(path) # delte zipped files
    print("Done!")

main(20, 'C:/Users/Omar2/Desktop/')

