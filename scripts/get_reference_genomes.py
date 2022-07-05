#!/usr/bin/env python

from ast import While
from fileinput import filename
from ftplib import FTP
from pickle import TRUE
import random
import os
import argparse
import os.path
import gzip
import shutil
from numpy import save

# login to ftp and direct to needed directory
try:
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login(user = "#####", passwd = "#####") # login creds to the ftp server
    ftp.cwd('genomes/.vol2/refseq/bacteria/') 
    ftp.encoding = "utf-8"
except:
    print("ALERT: Connection Failed! (check network or ftp credentials")

def goToDirect(): # go to the directory with the FASTA files
    final_Directory = False
    while final_Directory != True:
        for file in ftp.nlst():
            if file.endswith(".fna.gz"):
                final_Directory = True
        if final_Directory != True:
            ftp.cwd(ftp.nlst()[0])
    return final_Directory

def countFilesInDirectory(path): # count how many directories are in the reference_genemoes folder
    count = 0
    for directory in os.listdir(path):
        count += 1
    return count

def checkCount(path, INTIAL_COUNT, number_of_genemoes):
    count = 0
    for Directory in os.listdir(path):
        count += 1
    if count != INTIAL_COUNT + number_of_genemoes :
        return False
    return True

def downloadFNAfile(path, current_directory_name): # download all FASTA files in a directory
    path = makeFile(path + f"/{current_directory_name}")
    for filename in ftp.nlst():
        if filename.endswith(".fna.gz"):
            completeName = os.path.join(path, filename)         
            file = open(completeName, "wb")
            ftp.retrbinary(f"RETR {filename}", file.write)
            file.close()

def returnToOriginalDirect(): # return to main directory
    ftp.cwd('/') 
    ftp.cwd('genomes/.vol2/refseq/bacteria/') 

def makeFile(save_path): # make folders
    try: 
        os.mkdir(save_path) 
    except OSError as error: 
        print("ALERT: " + save_path + " directory currently exists, so the files will be downloaded or overlapped there!") 
    return save_path

def unZip(path): # unzip the files
    for bacteriaDirectory in os.listdir(path):
        new_path = path + f"/{bacteriaDirectory}"
        for fnaFILE in os.listdir(new_path):
            final_path = new_path + f"/{fnaFILE}"
            name, extension = os.path.splitext(fnaFILE)
            if extension == ".gz":
                try:
                    with gzip.open(final_path, 'rb') as f_in:
                        with open(final_path[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                except:
                    print(f"ALERT: The following file is too large and didnt unzip! \n*>* {final_path}")

def removeZip(path): # remove the old zip files
    for bacteriaDirectory in os.listdir(path):
        new_path = path + f"/{bacteriaDirectory}"
        for fnaFILE in os.listdir(new_path):
            name, extension = os.path.splitext(fnaFILE)
            if extension == ".gz":
                os.remove(os.path.join(new_path, fnaFILE))

def main():
    parser = argparse.ArgumentParser(description="This Script will connect to the ftp server,"
                                    " then will download a random sample of refrenced genemoes,"
                                    " you can run the script again to add more genemoes, "
                                    "to proceed you need to modify to your ftp crednetials, then"
                                    " enter the 3 postional arguments 'number_of_genemoes', "
                                    "'wanna_unzip', and 'save_path'"
                                    ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("number_of_genemoes", type=int, help="The number of refrence genemoes to be downloaded")
    parser.add_argument("wanna_unzip", type=str, help="Whether you would like the files to be unzipped after download."
                                                    "If files are too large, the unzipping operation fails", default= "True")
    parser.add_argument("save_path",type=str, help="The full path to the directory that the files will be downloaded in")
    args = parser.parse_args()
    number_of_genemoes = args.number_of_genemoes
    wanna_unzip = args.wanna_unzip
    save_path = args.save_path
    # ----------------------------- ^^^ Arguments and Interface ^^^ ---------------------------------------------
    print("Connecting...")
    bacteriaDirectoryNames = ftp.nlst() # list of possible refrence bacterias
    checker_varible = bacteriaDirectoryNames[0]
    randomlist = random.sample(range(0, len(bacteriaDirectoryNames)-2), number_of_genemoes) # get a random sample
    print("Downloading...")
    path = makeFile(save_path + 'refrence_genemoes') # make main folder
    INITAL_COUNT = countFilesInDirectory(path)
    used_genemoes = []
    while checkCount(path, INITAL_COUNT, number_of_genemoes) != True: # main loop to save file
        aRandomNumber = random.randint(0, len(bacteriaDirectoryNames)-2)
        while aRandomNumber in used_genemoes:
            aRandomNumber = random.randint(0, len(bacteriaDirectoryNames)-2)
        used_genemoes.append(aRandomNumber)
        try:
            try:
                current_directory_name = bacteriaDirectoryNames[aRandomNumber] # whats the bacteria name
                ftp.cwd(current_directory_name) # go to the specfied bacteria directory
                goToDirect() # go to the directory with the FASTA files
                downloadFNAfile(path, current_directory_name) # download the FASTA files
                returnToOriginalDirect() # return to the main directory once done
            except:  
                if ftp.nlst()[0] != checker_varible:
                    returnToOriginalDirect() # return to the main directory once done
        except EOFError as e:
            print(f"ALERT: Unkown Error occured! Please re-run for more genemoes if needed! \nTotal Number of downloaded genemoes: {countFilesInDirectory(path)}")
            break
    ftp.close() # close the connection
    if wanna_unzip == "true" or wanna_unzip == "True" or wanna_unzip == "yes" or wanna_unzip == "Yes":
        print("Unzipping files...")
        unZip(path) # unzip the files
        removeZip(path) # delte zipped files
    print("Done!")

if __name__ == "__main__":
    # issue1: explore why EOFError occurs often from ftplib moudle
    # issue2: enhance the method that goes to the propper directory of the geneome
    main()

