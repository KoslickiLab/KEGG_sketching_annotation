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

class Helper:
    def __init__(self, credentials):
        # Connect to the FTP database
        # Connect to FTP database
        # login to ftp and direct to needed directory
        try:
            ftp = FTP('ftp.ncbi.nlm.nih.gov')
            ftp.login(user="anonymous", passwd=f"{credentials}")  # login creds to the ftp server
            ftp.cwd('genomes/.vol2/refseq/bacteria/')
            ftp.encoding = "utf-8"
            self.ftp = ftp
        except:
            raise Exception(
                "ALERT: Could not connect to the ftp server, please check your internet connection or ftp credentials")

    def go_to_direct(self): # go to the directory with the FASTA files
        ftp = self.ftp
        final_Directory = False
        while final_Directory != True:
            for file in ftp.nlst():
                if file.endswith(".fna.gz"):
                    final_Directory = True
            if final_Directory != True:
                ftp.cwd(ftp.nlst()[0])
        return final_Directory

    def count_files_in_directory(self, path): # count how many directories are in the reference_genomes folder
        count = 0
        for directory in os.listdir(path):
            count += 1
        return count

    def check_count(self, path, INTIAL_COUNT, number_of_genomes):
        count = 0
        for Directory in os.listdir(path):
            count += 1
        if count != INTIAL_COUNT + number_of_genomes:
            return False
        return True

    def download_FNA_file(self, path, current_directory_name): # download all FASTA files in a directory
        ftp = self.ftp
        path = self.make_file(os.path.join(path, f"{current_directory_name}"))
        for filename in ftp.nlst():
            if filename.endswith(".fna.gz"):
                completeName = os.path.join(path, filename)
                file = open(completeName, "wb")
                ftp.retrbinary(f"RETR {filename}", file.write)
                file.close()

    def return_to_original_direct(self): # return to main directory
        ftp = self.ftp
        ftp.cwd('/')
        ftp.cwd('genomes/.vol2/refseq/bacteria/')

    def make_file(self, save_path): # make folders
        try:
            os.mkdir(save_path)
        except OSError as error:
            print("ALERT: " + save_path + " directory currently exists, so the files will be downloaded or overlapped there!")
        return save_path

    def unzip(self, path):  # unzip the files
        for bacteriaDirectory in os.listdir(path):
            new_path = os.path.join(path, f"{bacteriaDirectory}")
            for fnaFILE in os.listdir(new_path):
                final_path = os.path.join(new_path, f"{fnaFILE}")
                name, extension = os.path.splitext(fnaFILE)
                if extension == ".gz":
                    try:
                        with gzip.open(final_path, 'rb') as f_in:
                            with open(final_path[:-3], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    except:
                        print(f"ALERT: The following file is too large and didnt unzip! \n*>* {final_path}")

    def remove_zip(self, path): # remove the old zip files
        for bacteriaDirectory in os.listdir(path):
            new_path = os.path.join(path, f"{bacteriaDirectory}")
            for fnaFILE in os.listdir(new_path):
                name, extension = os.path.splitext(fnaFILE)
                if extension == ".gz":
                    os.remove(os.path.join(new_path, fnaFILE))

def main():
    parser = argparse.ArgumentParser(description="This Script will connect to the ftp server,"
                                    " then will download a random sample of refrenced genomes,"
                                    " you can run the script again to add more genomes, "
                                    "to proceed you need to modify to your ftp crednetials, then"
                                    " enter the 3 postional arguments 'number_of_genomes', "
                                    "'wanna_unzip', and 'save_path'"
                                    ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-n", "--number_of_genomes", type=int, help="The number of refrence genomes to be downloaded")
    parser.add_argument("-u", "--unzip", action='store_true', help="Whether you would like the files to be unzipped after download."
                                                    "If files are too large, the unzipping operation fails")
    parser.add_argument("-s", "--save_path", type=str, help="The full path to the directory that the files will be downloaded in")
    parser.add_argument("-c", "--credentials", type=str, help="NCBI wants you to login to the ftp server with your email address", default="user@email")
    args = parser.parse_args()
    number_of_genomes = args.number_of_genomes
    wanna_unzip = args.unzip
    save_path = args.save_path
    credentials = args.credentials
    # Check arguments
    if number_of_genomes <= 0:
        raise Exception("The number of genomes to be downloaded must be greater than 0")
    if not os.path.exists(save_path):
        print(f"The file path {save_path} does not exist, creating it now")
        os.makedirs(save_path)
    if not os.path.isdir(save_path):
        raise Exception(f"The file path must be a directory. I was given {save_path}")
    # check if directory is empty
    if os.listdir(save_path):
        print(f" WARNING: the directory {save_path} is not empty, press any key to continue or ctrl+c to exit")
        input()

    # Instantiate the helper class
    helper = Helper(credentials)

    # ----------------------------- ^^^ Arguments and Interface ^^^ ---------------------------------------------
    print("Connecting...")
    # Keep trying to connect until we succeed
    got_names = False
    bacteriaDirectoryNames = None
    while not got_names and not bacteriaDirectoryNames:
        try:
            bacteriaDirectoryNames = helper.ftp.nlst()  # list of possible refrence bacterias
            got_names = True
        except EOFError:
            print("FTP connection error, trying again")
    checker_varible = bacteriaDirectoryNames[0]
    randomlist = random.sample(range(0, len(bacteriaDirectoryNames)-2), number_of_genomes)  # get a random sample
    print("Downloading...")
    path = helper.make_file(save_path + 'refrence_genomes')  # make main folder
    INITAL_COUNT = helper.count_files_in_directory(path)
    used_genomes = []
    while helper.check_count(path, INITAL_COUNT, number_of_genomes) != True:  # main loop to save file
        aRandomNumber = random.randint(0, len(bacteriaDirectoryNames)-2)
        while aRandomNumber in used_genomes:
            aRandomNumber = random.randint(0, len(bacteriaDirectoryNames)-2)
        used_genomes.append(aRandomNumber)
        try:
            try:
                current_directory_name = bacteriaDirectoryNames[aRandomNumber]  # whats' the bacteria name
                helper.ftp.cwd(current_directory_name)  # go to the specified bacteria directory
                helper.go_to_direct()  # go to the directory with the FASTA files
                helper.download_FNA_file(path, current_directory_name)  # download the FASTA files
                helper.return_to_original_direct()  # return to the main directory once done
            except:  
                if helper.ftp.nlst()[0] != checker_varible:
                    helper.return_to_original_direct()  # return to the main directory once done
        except EOFError as e:
            print(f"ALERT: Unkown Error occured! Please re-run for more genomes if needed! \nTotal Number of downloaded genomes: {helper.count_files_in_directory(path)}")
            break
    helper.ftp.close() # close the connection
    if wanna_unzip == "true" or wanna_unzip == "True" or wanna_unzip == "yes" or wanna_unzip == "Yes":
        print("Unzipping files...")
        helper.unzip(path)  # unzip the files
        helper.remove_zip(path)  # delete zipped files
    print("Done!")

if __name__ == "__main__":
    # issue1: explore why EOFError occurs often from ftplib moudle
    # issue2: enhance the method that goes to the propper directory of the geneome
    main()

