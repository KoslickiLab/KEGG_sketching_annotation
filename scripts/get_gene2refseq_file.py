#!/usr/bin/env python

import os
import argparse
from ftplib import FTP

def main():
    parser = argparse.ArgumentParser(description="This Script will connect to the ftp server,"
                                    " then will download the gene2refseqfile, to procceed you need to"
                                    " modify the username and password of the ftp login (line 16)."
                                    ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-u', '--user_name', type=str, help="User name for the ftp login. NCBI wants an email address.",
                        default="user@email.com")
    args = parser.parse_args()
    user_name = args.user_name
    # login to ftp and direct to needed directory
    print("Connecting...")
    try:
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login(user="anonymous", passwd=f"{user_name}")  # login creds to the ftp server
        ftp.cwd('gene/DATA/') 
        ftp.encoding = "utf-8"
        filename = "gene2refseq.gz"
    except:
        raise Exception("ALERT: Connection Failed! (check network or ftp credentials")
    print("downloading...")
    try:
        ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)  # download needed file
    except:
        raise Exception("ALERT: Download Failed! (check disk space or connection)")
    ftp.quit()  # quit connection
    print("Done!")


if __name__ == "__main__":
    main()
