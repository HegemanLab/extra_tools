'''
Little tool that takes in list of m/zs and then outputs
each mz that was present as well as how many times it showed up
overall in all of the inputted lists of m/zs. Rounds to 4 decimal places
as is.
'''


import csv
import pandas as pd


# short helper to pass into map. Takes number and rounds to 3 decimal places
def roundto4(x):
    return round(x, 4)

# Takes list of values and returns the unique m/z values and their counts
def finder(in_list):

    # Rounds each value in each list of mzs (1 list per file input)
    for i in range(len(in_list)):
         in_list[i] = map(roundto4, in_list[i])

    # Removes duplicates
    no_dups = []
    for i in in_list:
        for j in i:
            if j not in no_dups:
                no_dups.append(j)

    counts = [0] * len(no_dups)

    # For each unique m/z, get a count of how many times it appears
    for x in no_dups:
        for i in in_list:
            if x in i:
                counts[no_dups.index(x)] += 1

    return [no_dups, counts]


# Takes a list containing two lists of equal lengths and writes them to
# a csv (need to have .csv in filename
def writeToCSV(mzs_and_counts, output_filename):
    with open(output_filename, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['mz', 'count'])
        for x in range(0, len(mzs_and_counts[1])):
            writer.writerow([mzs_and_counts[0][x], mzs_and_counts[1][x]])


# Reads in mz values from a file and returns list of mzs
# FILE MUST HAVE 'mz' AS COLUMN HEADER
def readFromCSV(filename):
    df = pd.read_csv(filename)
    return list(df['mz']) # Make sure file has a column named mz


def processDuplicates(files_list, outputname):
    mzs = []
    for f in files_list:
        mzs.append(readFromCSV(f))
    writeToCSV(finder(mzs), outputname)


# # Example use
# nfiles = [
#     'neg-ACM_sept16_T1R2_GL2_method1.csv',
#     'neg-ACM_sept16_T1R2_GL7_method1.csv',
#     'neg-ACM_sept16_T1R3_GL7_method1.csv',
#     'neg-ACM_sept16_T1R3_GL21_method1.csv'
# ]
#
# pfiles = [
#     'pos-ACM_sept16_T1R2_GL2_method1.csv',
#     'pos-ACM_sept16_T1R2_GL7_method1.csv',
#     'pos-ACM_sept16_T1R3_GL7_method1.csv',
#     'pos-ACM_sept16_T1R3_GL21_method1.csv'
# ]
#
# processDuplicates(nfiles, 'neg-ACM-dups.csv')
# processDuplicates(pfiles, 'pos-ACM-dups.csv')

print 'boom'

