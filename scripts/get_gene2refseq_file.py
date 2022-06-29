import os
from ftplib import FTP

def main():
    # login to ftp and direct to needed directory
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login(user = "anonymous", passwd = "ohr5003") # login creds to the ftp server
    ftp.cwd('gene/DATA/') 
    ftp.encoding = "utf-8"
    filename = "gene2refseq.gz"
    print("downloading...")
    ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)
    ftp.quit()
    print("done !")

main()