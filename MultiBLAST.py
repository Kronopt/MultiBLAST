#!python2
# coding: utf-8

"""
MultiBLAST
Sequentially Blast every fasta file (.faa or .fas) on a folder against a local database

DEPENDENCIES:
    - Python 2.7
    - Blast+ executables (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - LINUX (tested and developed on Ubuntu 14.0.4, but will most likely work on WINDOWS and MAC too)

HOW TO RUN:
    - Organize your files as follows on a folder:
        - MultiBLAST.py
        - database files (obtained via makeblast, with 'makeblastdb -in [yourGenes].pfasta -dbtype prot')
        - folder with fasta files (.faa or .fas) to blast (must be the only existing folder in the directory)
    - Run MultiBLAST.py from the command-line as follows 'python MultiBLAST.py [-t [#]]'
    - For more help, type 'python MultiBLAST.py -h'

ENSURES:
    - Each fasta file will be blasted against the given database
    - Each blast output is created inside the fasta folder
"""

import os
import argparse
from time import time, strftime, localtime
from sys import exit

__author__ = 'Pedro HC David, https://github.com/Kronopt'
__credits__ = ['Pedro HC David']
__version__ = '1.0'
__date__ = '00:17h, 21/10/2016'
__status__ = 'Finished'


def blast(threads):
    """
    Main Function

    PARAMETERS:
        threads : int/str
            Number of CPU threads to use

    REQUIRES:
        threads must be a number >= 1

    ENSURES:
        Blast each file against the given database
    """

    root_folder, root_list = root_files()

    # Obtains root path, database name and a list containing all the fasta failes
    folder_root_dir = fasta_folder(root_list)
    db = db_file(root_list)
    fasta_files = fasta_files_folder(root_folder, folder_root_dir)

    start_blast_text()
    
    # Timer, START
    start = time()

    # Time of program start
    program_start_time = strftime('%d/%m/%y  %H:%M:%S', localtime())

    # Main loop
    loop_blast(fasta_files, root_folder, folder_root_dir, db, threads)
    end_blast_text(start, program_start_time)


def exit_script():
    """
    Exits the script
    """
    
    print
    print '*****'
    print 'EXITING SCRIPT'
    print '*****'
    exit()


# -----------------------------------------------------
# FILES AND FOLDERS SECTION
# -----------------------------------------------------


def root_files():
    """
    Creates a list of every file/folder in the current working directory (ROOT)
    """
    
    root_folder = os.getcwd()
    root_list = os.listdir(root_folder)
    root_list.sort()

    return root_folder, root_list


def fasta_folder(root_list):
    """
    Name of the folder containing the fasta files (.faa or .fas) files
    Prints the folder to the screen

    PARAMETERS:
        root_list : list of str
            List of directories in the root directory

    REQUIRES:
        Directories must be valid

    ENSURES:
        Path to the only folder in the root directory
    """

    print
    
    folder_root_dir = filter(lambda x: os.path.isdir(x), root_list)
    if len(folder_root_dir) != 1:
        print 'ERROR: None or more than one folder was found'
        exit_script()

    print 'FASTA FOLDER:'
    print '    ', folder_root_dir[0]
    print

    return folder_root_dir[0]


def db_file(root_list):
    """
    Database files
    Prints the name of the database to the screen

    PARAMETERS:
        root_list : list of str
            List of directories in the root directory

    REQUIRES:
        Directories must be valid

    ENSURES:
        Path to the .pfasta database file
    """
    
    db = filter(lambda x: x.endswith('.pfasta'), root_list)
    if len(db) != 1:
        print 'ERROR: None or more than one database file was found'
        exit_script()

    print 'DB FILE:'
    print '    ', db[0]
    print

    return db[0]


def fasta_files_folder(root_folder, folder_root_dir):
    """
    Fasta files
    Prints every file name to the screen

    PARAMETERS:
        root_folder : str
            Root folder's path
        folder_root_dir : str
            Fasta file folder's name

    REQUIRES:
        root_folder must be a valid directory
        folder_root_dir must be a valid folder name

    ENSURES:
        List of fasta files
    """

    fasta_files = filter(lambda x: (x.endswith('.faa') or x.endswith('.fas'))
                         and os.path.isfile(root_folder + '/' + folder_root_dir + '/' + x),
                         os.listdir(root_folder + '/' + folder_root_dir))
    fasta_files.sort()

    print 'FASTA FILES:'
    for faa_file in fasta_files:
        print '    ', faa_file
    print

    return fasta_files


# -----------------------------------------------------
# blast SECTION
# -----------------------------------------------------


def start_blast_text():
    """
    Start blast text
    """
    
    print
    print '*****'
    print 'STARTING BLASTS'
    print '*****'
    print
    print


def loop_blast(fasta_files, root_folder, folder_root_dir, db, threads):
    """
    Main code for blasting each fasta file
    Prints information on each successful blast

    PARAMETERS:
        fasta_files : list of str
            List of fasta files
        root_folder : str
            Path to root folder
        folder_root_dir : str
            Fasta folder's name
        db : str
            Database name
        threads : str/int
            Number of CPU threads to use

    REQUIRES:
        fasta_files must have a list of existing fasta files
        root_folder must be a valide directory
        folder_root_dir must be the name of an existing folder
        db must be the name of the existing database
        threads must be a number >= 1

    ENSURES:
        Blast of each fasta file agains the local specified database
    """

    count = 1
    for fasta_file in fasta_files:

        # Timer, START
        start = time()
        
        print '------------------'
        print 'blast', count, 'out of', len(fasta_files)
        print 'BLASTING:'
        print '    ', fasta_file

        out_file_name = fasta_file[:len(fasta_file) - 4] + '__' + db[0:-7] + '_OUT'

        # Tries to blast
        standard_blast = 'blastp -out ' + out_file_name + ' -outfmt 6 -query "' + root_folder + '/' + \
                         folder_root_dir + '/' + fasta_file + '" -db ' + db + ' -num_threads ' + str(threads)
        
        try:
            os.system(standard_blast)
        except:
            print
            print 'ERROR: Could not run NCBI blast+ program blastp'
            exit_script()

        # Timer, END
        end = time()

        days, hours, minutes, seconds = sec_to_hours(end - start)

        print
        print 'TIME ELAPSED:'
        if days != 0:
            print str(days) + 'd ',
        print str(hours) + 'h', str(minutes) + 'm', str(seconds) + 's'
        print

        # Moves the output file to the output file folder
        os.rename(out_file_name, folder_root_dir + '/' + out_file_name)

        count += 1


def end_blast_text(timer_start, program_start):
    """
    End blast text

    PARAMETERS:
        timer_start : float
            Value obtained from time.time()
        program_start : str
            Current time in text format

    ENSURES:
        Time information printed to the screen
    """

    days, hours, minutes, seconds = sec_to_hours(time() - timer_start)
    program_end = strftime('%d/%m/%y  %H:%M:%S', localtime())
    
    print '------------------'
    print 'TIME OF START/END:'
    print program_start
    print program_end
    print
    print 'DURATION:'
    if days != 0:
        print str(days) + 'd ',
    print str(hours) + 'h', str(minutes) + 'm', str(seconds) + 's'
    print
    print
    print '*****'
    print 'END OF BLASTS'
    print '*****'
    print


def sec_to_hours(seconds):
    """
    Converts seconds into hours (or days, if it takes that long to run)

    PARAMETERS:
        seconds : int/str/float
            number of seconds

    REQUIRES:
        seconds must be a number >= 1

    ENSURES:
        Seconds translated into days, hours, minutes and seconds
    """

    seconds = int(seconds)
    
    days = seconds / 86400
    hours = (seconds - days * 86400) / 3600
    minutes = (seconds - days * 86400 - hours * 3600) / 60
    seconds = (seconds - days * 86400 - hours * 3600 - minutes * 60)

    return days, hours, minutes, seconds


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Blasts fasta files (.faa or .fas) to a local database. '
                                                 'Organize your files so that the root folder has only: [1] this '
                                                 'script (MultiBLAST.py); [2] all the database files (obtained via '
                                                 '"makeblastdb -in [yourGenes].pfasta -dbtype prot"); [3] folder with '
                                                 'fasta files (.faa or .fas) to blast (must be the only folder '
                                                 'available)')
    parser.add_argument('-t', metavar='#', nargs='?', default=1, type=int, const=1, help='Number of CPU threads')
    arguments = parser.parse_args()

    blast(arguments.t)
