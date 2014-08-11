#------------------------------------------------------------------------------#
#                              popsites.py                                     #
#------------------------------------------------------------------------------#

# This script converts a list of SNP sites or genes .txt format to a new tabix
# formatted .txt file and feeds it into tabix using xargs. It includes the 
# option to add a buffer of a specified number of base pairs to either side
# of each site or gene.

# Running txt_to_tabix.py requires that both xargs and tabix be available.

# Usage: python3 pop_sites.py [-h] [-b BUFF] -f FDR [-g GTF] -i INPATH -p POP
#                             (-c | -r REP)


# Firstly, import dependencies.
import argparse
import csv
import ftplib
import os
import random
import numpy
from scipy import stats



#------------------------------------------------------------------------------#
#                                Arguments                                     #
#------------------------------------------------------------------------------#

# Initialize the argument parser and parse arguments.

parser = argparse.ArgumentParser( 
        description='''Identify population-specific SNP sites from a given 
                       input set of SNP sites or genes.''')
parser.add_argument('-b','--buff', required=False, type=int, default=0,
        help='Size of buffer around sites or genes.')
parser.add_argument('-f','--fdr', required=True, type=float,
        help='False discovery rate for population-specific SNP sites.')
parser.add_argument('-g','--gtf', required=False,
        help='Path to the .gtf file used for gene-coordinate conversion.')
parser.add_argument('-i','--inpath', required=True,
        help='Path to the input .txt containing slicing data.')
parser.add_argument('-p','--pop', required=True,
        help='Population or superpopulation to test, e.g. CEU or EUR.')
testMEG = parser.add_mutually_exclusive_group(required=True)
testMEG.add_argument('-c','--chisq', action='store_true',
        help='Use a chi-squared test to identify population-specific SNPs.')
testMEG.add_argument('-r','--rep', type=int,
        help='''Number of replicates to permute when constructing empirical
                probability distributions.''')
args = parser.parse_args()




#------------------------------------------------------------------------------#
#                          Functional Definitions                              #
#------------------------------------------------------------------------------#

# Here we define functions to carry out the following tasks:
# 1. Subsetting the population panel.
# 2. Parsing site inputs and converting them to tabix input format.
# 3. parsing gene inputs and converting them to tabix input format.
# 4. Executing tabix on the formatted inputs.
# 5. Subsetting the population-specific VCF data.
# 6. Computing allele frequencies with vcftools.
# 7. Identifying population-specific SNPs and genes via chi-squared test.
# 8. Constructing empirical frequency distributions.
# 9. Identifying population-specific SNPs and genes via permutation test.
# 10. Cleaning up.




#------------------------------------------------------------------------------#
#                     1. Subsetting the population panel                       #
#------------------------------------------------------------------------------#

# Generate a new panel file including only data from the indicated population.

def subset_pop(pop):
    print('Subsetting the population panel')

    # Make sure the input population is valid
    
    if pop not in '''AFR AMR ASN EUR SAN MSL EXN ASW ACB MXL PUR CLM PEL CHB JPT
                     CHS CDX CEU TSI FIN GBR IBS PJL BEB STU ITU''':
        raise Exception('Not a valid input population')
    
    # Indicate the locations where the full and subsetted population panels
    # will be stored.
    
    all_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', 'ALL', '.panel')
    pop_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', pop, '.panel')
    
    # Parse superpopulation inputs.
    
    pop = pop.replace('AFR', 'MSL ESN ASW ACB'
            ).replace('AMR', 'MXL PUR CLM PEL'
            ).replace('ASN', 'CHB JPT CHS CDX'
            ).replace('EUR', 'CEU TSI FIN GBR IBS'
            ).replace('SAN', 'PJL BEB STU ITU'
            )
            
    # Connect to the 1000 genomes EBI FTP site and retrieve the population
    # panel.

    ftp = ftplib.FTP('ftp.1000genomes.ebi.ac.uk')
    ftp.login()
    f = open(all_panel, 'wb')
    ftp.retrbinary('RETR /vol1/ftp/release/20100804/20100804.ALL.panel',
                  f.write
                 )
    f.close()
    ftp.quit()
    
    # Initialize a file to store the input population subject ID's
    # and a variable to check for empty intersection between the input set of
    # SNPs and the set detected in the input population.
    
    empty_intersection = True
    output = open(pop_panel, 'w')
    
    # Read the retrieved panel file line by line and store the lines referring
    # to the EUR samples.
    
    with open(all_panel) as f:
        for line in f:
            parsed_line = line.split('\t')
            if parsed_line[1] in pop:
                output.write(line)
                empty_intersection = False
                
    # Close the output file.
    
    output.close()
    
    # Check for empty intersection
    
    if empty_intersection:
        os.system(' '.join(['rm', all_panel, pop_panel]))
        raise Exception('Input sites not detected in input population.')
        
        
        
        
#------------------------------------------------------------------------------#
#       2. Parsing site inputs and converting them to tabix input format.      #
#------------------------------------------------------------------------------#

def convert_sites(inpath, buff):
    print('Converting sites to tabix input format')

    # Initialize the output file.
    
    output = open(inpath[::-1][4:][::-1]+'_tabix.txt', 'w')
    
    # Read the input file line by line and record an interval centered on each
    # indicated site in tabix input format.
    
    with open(inpath) as f:
        for line in f:
            parsed_line = line.split('\t')[0:2]
            try:
                output.write('{}:{}-{}{}'.format(parsed_line[0],
                                             str(int(parsed_line[1])-buff),
                                             str(int(parsed_line[1])+buff+1),
                                             '\n'
                                            ))
            except ValueError:
                pass
                
   # Close the output file.

    output.close()
    
    
    
    
#------------------------------------------------------------------------------#
#       3. Parsing gene inputs and converting them to tabix input format.      #
#------------------------------------------------------------------------------#

def convert_genes(inpath, gtf, buff):
    print('Converting genes to tabix input format')
    
    # Collect a dictionary keying each gene to its absolute coordinates as 
    # retrieved from the gtf file, with a specified buffer on either side.

    hg19 = {}

    with open(gtf) as f:
        for line in f:
            try:
                parsed_line = line.split('\t')
                chrom = parsed_line[0].replace('chr','')
                start = str(int(parsed_line[3])-buff)
                end = str(int(parsed_line[4])+buff+1)
                gene = parsed_line[8].split(' ')[1].replace('"','').replace(';','')
                try:
                    hg19[gene] = [chrom,
                                  str(min(int(hg19[gene][1]),int(start))),
                                  str(max(int(hg19[gene][2]),int(end)))
                                 ]
                except KeyError:
                    hg19[gene] = [chrom, start, end]
            except IndexError:
                pass
            
    # Initialize the output file.

    output = open(inpath[::-1][4:][::-1]+'_tabix.txt', 'w')

    # Read the gtf file line by line to collect a dictionary keying gene ID's to
    # their absolute coordinates with a 100-bp cushion on either side.

    with open(inpath) as f:
        for line in f:
            gene = line.split('\t')[1].replace('\n','')
            try:
                output.write('{}:{}-{}{}'.format(hg19[gene][0],
                                                 hg19[gene][1],
                                                 hg19[gene][2],
                                                 '\n'
                                                ))
            except KeyError:
                pass
    
    # Close the output file.

    output.close()




#------------------------------------------------------------------------------#
#                4. Executing tabix on the formatted inputs.                   #
#------------------------------------------------------------------------------#

def execute_tabix(inpath):
    print('Executing tabix')
    
    # Identify the .vcf file to be used.
    
    vcf = ''.join(['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/',
                   'ALL.2of4intersection.20100804.genotypes.vcf.gz'
                  ])
    
    # Use tabix to generate a header.
    
    os.system(' '.join(['tabix',
                        '-fh', vcf, 'chr1:1','>', inpath[::-1][4:][::-1]+'.vcf'
                       ]))
    
    # Call tabix on the converted list of genomic regions, using xargs to 
    # feed in the inputs.
                          
    os.system(' '.join(['xargs',
                        '-a', inpath[::-1][4:][::-1]+'_tabix.txt','-I', '{}',
                        'tabix',
                        '-f', vcf, '{}','>>', inpath[::-1][4:][::-1]+'.vcf'
                       ]))

    
    
    
#------------------------------------------------------------------------------#
#              5. Subsetting the population-specific VCF data.                 #
#------------------------------------------------------------------------------#

def subset_vcf(inpath, pop):

    print('Subsetting the .vcf file with vcf-subset.')

    # Indicate the location of the full and subsetted population panels.
    
    all_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', 'ALL', '.panel')
    pop_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', pop, '.panel')
    
    # Initialize a data structure to store a list of samples from the input
    # population.
    
    all_samples = []
    pop_samples = []
    
    # Read the both panels line by line and record their respective lists of
    # samples.
    
    with open(all_panel) as f:
        for line in f:
            all_samples.append(line.split('\t')[0])
    all_samples = ','.join(all_samples)
    
    with open(pop_panel) as f:
        for line in f:
            pop_samples.append(line.split('\t')[0])
    pop_samples = ','.join(pop_samples)
    
    # Feed the list of samples to vcf-subset, generating a full .vcf file with
    # NAN frequencies ommitted and a subsetted .vcf file containing data on the
    # input population only.
    
    os.system(' '.join(['vcf-subset',
                        '-c', all_samples,
                        inpath[::-1][4:][::-1]+'.vcf',
                        '>', inpath[::-1][4:][::-1]+'_'+'ALL'+'.vcf'
                       ]))
    os.system(' '.join(['vcf-subset',
                        '-c', pop_samples,
                        inpath[::-1][4:][::-1]+'.vcf',
                        '>', inpath[::-1][4:][::-1]+'_'+pop+'.vcf'
                       ]))
               
               
               
                       
#------------------------------------------------------------------------------#
#               6. Computing allele frequencies with vcftools.                 #
#------------------------------------------------------------------------------#

def compute_frequencies(inpath, pop):
    print('Computing allele frequencies with vcftools')
    
    # Sort the full .vcf file and compute allele frequencies.
    
    os.system(' '.join(['vcf-sort', 
                        inpath[::-1][4:][::-1]+'_'+'ALL'+'.vcf',
                        '>', 
                        inpath[::-1][4:][::-1]+'_'+'ALL'+'_sort.vcf'
                       ]))
    os.system(' '.join(['vcftools',
                        '--vcf', inpath[::-1][4:][::-1]+'_'+'ALL'+'_sort.vcf',
                        '--freq',
                        '--out', inpath[::-1][4:][::-1]+'_'+'ALL'
                       ]))
                       
    # Sort the population-subsetted .vcf file and compute allele frequencies.                   
                       
    os.system(' '.join(['vcf-sort', 
                        inpath[::-1][4:][::-1]+'_'+pop+'.vcf',
                        '>', 
                        inpath[::-1][4:][::-1]+'_'+pop+'_sort.vcf'
                       ]))
    os.system(' '.join(['vcftools',
                        '--vcf', inpath[::-1][4:][::-1]+'_'+pop+'_sort.vcf',
                        '--freq',
                        '--out', inpath[::-1][4:][::-1]+'_'+pop
                       ]))




#------------------------------------------------------------------------------#
#   7. Identifying population-specific SNPs and genes via chi-squared test.    #
#------------------------------------------------------------------------------#

def chi_squared_test(fdr, inpath, pop):
    print('Applying parametric chi-squared test')
    
    # Identify a few important filenames.
    
    all_frq = inpath[::-1][4:][::-1]+'_'+'ALL'+'.frq'
    pop_frq = inpath[::-1][4:][::-1]+'_'+pop+'.frq'
    dist_path = inpath[::-1][4:][::-1]+'_dist.csv'
    verdict_path = inpath[::-1][4:][::-1]+'_verdict.txt'
            
    # Construct a dictionary keying each site listed in all_frq (grand
    # population frequency file) to its computed frequency values.
    
    all_dict = {}
    all_dict['CHROM:POS'] = '{ALLELE:FREQ}'
    with open(all_frq) as f:
        for line in f:
            parsed_line = line.split('\t')
            try:
                all_dict[':'.join([parsed_line[0],
                                   parsed_line[1]
                                  ])] = [float(j[2:]) for j in [
                                                 parsed_line[4],
                                                 parsed_line[5].replace('\n','')
                                                               ]]
            except IndexError:
                pass
    
    # Initialize a table to store sites  and a counter for the number of
    # sites listed in pop_frq.
    
    p_vals = []
    pop_sites = []
    n_sites = 0
    
    # For each of the pop_frq sites, compute a chi-squared statistic and record
    # a p-value. Also, count the number of sites listed in pop_frq and construct
    # a list of them.
    
    with open(pop_frq) as f:
        next(f)
        for line in f:
            n_sites += 1
            parsed_line = line.split('\t')
            site = ':'.join([parsed_line[0], parsed_line[1]])
            
            # Compute chi-squared statistic and p-value
            
            observed = numpy.array([float(j[2:]) 
                                    for j in 
                                    [parsed_line[4],
                                     parsed_line[5].replace('\n','')
                                    ]])
            expected = numpy.array(all_dict[site])
            p_val = stats.chisquare(observed, expected)[1]
            
            # Append this result to the p_vals table
        
            p_vals.append([site, p_val])
            
    # Carry out a Benjamini-Hochberg procedure on the recorded p-values and
    # tag each site as 'POP-SPECIFIC' or 'INCONCLUSIVE'. Then write the results
    # to a tab-delimited output file.
    
    # Sort the sites by p-value
    
    p_vals.sort( key=lambda entry: entry[1])
    
    # Apply Benjamini-Hochberg
    
    cutoff = 0
    for i in list(range(n_sites)):
        if p_vals[i][1] <= i/n_sites * fdr:
            cutoff = int(i)
    for i in list(range(n_sites)):
        if i <= cutoff:
            p_vals[i].append(pop+'-SPECIFIC')
        else:
            p_vals[i].append('INCONCLUSIVE')
            
    # Sort the sites by coordinates.        
    
    p_vals.sort( key=lambda entry: entry[0])
    
    # Write the results to output.
    
    with open(verdict_path, 'w') as f:
        f.write('\t'.join(['SITE', 'P_VALUE','VERDICT']))
        for entry in p_vals:
            f.write('\t'.join([str(x) for x in entry])+'\n')
            
    # Print a summary of the results.
    
    print(''.join(['Analysis complete: discovered ', str(cutoff), ' ', pop, 
                   '-specific sites and ', str(n_sites-cutoff),
                   ' inconclusive sites.'
                  ]))
    



#------------------------------------------------------------------------------#
#             8. Constructing empirical frequency distributions.               #
#------------------------------------------------------------------------------#

def empirical_dist(inpath, pop, rep):
    print('Sampling for empirical frequency distributions.')
    
    # Identify a few important filenames.
    
    all_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', 'ALL', '.panel')
    pop_panel = '{}/{}{}{}'.format(os.getcwd(), '20100804.', pop, '.panel')
    pop_frq = inpath[::-1][4:][::-1]+'_'+pop+'.frq'
    sorted_vcf = inpath[::-1][4:][::-1]+'_ALL_sort.vcf'
    random_vcf = inpath[::-1][4:][::-1]+'_rand.vcf'
    random_frq = inpath[::-1][4:][::-1]+'_rand.frq'
    dist_path = inpath[::-1][4:][::-1]+'_dist.csv'
    
    # Collect a list of all of the 1000 genomes samples, which will be used as
    # a pool from which to draw a random subset.
            
    all_samples = []
    with open(all_panel) as f:
        for line in f:
            all_samples.append(line.split('\t')[0])
            
    # Count the number of samples in the input population, so that the random
    # subsets drawn from the grand population will be of equal size.
            
    pop_size = 0
    with open(pop_panel) as f:
        for line in f:
            pop_size +=1
            
    # Collect a list of all sites in the pop_frq file to be used as a header
    # for the empirical distribution file.
            
    pop_sites = []
    with open(pop_frq) as f:
        for line in f:
            parsed_line = line.split('\t')
            pop_sites.append(':'.join([parsed_line[0], parsed_line[1]]))
            
            
   # Initialize the empirical distribution file and write the header.

    with open(dist_path, 'w', newline='') as csv_file:
        dist_writer = csv.writer(csv_file, 
                                delimiter = ',',
                                quotechar='"',
                                quoting=csv.QUOTE_MINIMAL
                               )
        dist_writer.writerow(pop_sites)
        
        # A number of times equal to the value given for rep...
        
        for i in range(rep):
        
            # Choose a random subset of the grand population that is equal in 
            # size to the input population.
        
            random_subset = ','.join(random.sample(all_samples, pop_size))
            
            # Extract the randomly chosen subset from the full .vcf file with 
            # vcf-subset.
            
            os.system(' '.join(['vcf-subset', 
                                '-c', random_subset,
                                sorted_vcf,
                                '>', random_vcf
                               ]))
                               
            # Use vcftools to compute allele frequencies on the randomly 
            # subsetted .vcf file.
            
            os.system(' '.join(['vcftools', 
                                '--vcf', random_vcf, 
                                '--freq',
                                '--out', inpath[::-1][4:][::-1]+'_rand'
                               ]))
            
            # The set of sites for which frequencies are computed in the 
            # random_frq file may not be the same as the set computed in the 
            # pop_frq file. This can cause problems when constructing the 
            # empirical distribution file. To account for it, we collect a
            # dictionary keying each site listed in the random_frq file to its
            # computed freqency data.
            
            random_dict = {}
            random_dict['CHROM:POS'] = '{ALLELE:FREQ}'
            with open(random_frq) as f:
                for line in f:
                    parsed_line = line.split('\t')
                    try:
                        random_dict[':'.join([parsed_line[0],
                                              parsed_line[1]
                                             ])] = ' '.join([parsed_line[4],
                                                             parsed_line[5]
                                                            ]).replace('\n','')
                    except IndexError:
                        pass
            
            # Write a row of the empirical distribution file using the
            # dictionary collected in the previous step. Enter "." for missing
            # values.
            
            dist_row = []
            for key in pop_sites:
                try:
                    dist_row.append(random_dict[key])
                except KeyError:
                    dist_row.append('.')
            dist_writer.writerow(dist_row)            




#------------------------------------------------------------------------------#
#   9. Identifying population-specific SNPs and genes via permutation test.    #
#------------------------------------------------------------------------------#

def permutation_test(fdr, inpath, pop, rep):
    print('Carrying out permutation tests.')
    
    # Identify a few important filenames.
    
    all_frq = inpath[::-1][4:][::-1]+'_'+'ALL'+'.frq'
    pop_frq = inpath[::-1][4:][::-1]+'_'+pop+'.frq'
    dist_path = inpath[::-1][4:][::-1]+'_dist.csv'
    verdict_path = inpath[::-1][4:][::-1]+'_verdict.txt'
            
    # Construct a dictionary keying each site listed in all_frq (grand
    # population frequency file) to its computed frequency values.
    
    all_dict = {}
    all_dict['CHROM:POS'] = '{ALLELE:FREQ}'
    with open(all_frq) as f:
        for line in f:
            parsed_line = line.split('\t')
            try:
                all_dict[':'.join([parsed_line[0],
                                   parsed_line[1]
                                  ])] = [float(j[2:]) for j in [
                                                 parsed_line[4],
                                                 parsed_line[5].replace('\n','')
                                                               ]]
            except IndexError:
                pass
    
    # Construct a dictionary keying each site listed in pop_frq to its computed
    # frequency values. Also collect a list of sites listed in pop_frq and count
    # them.
    
    pop_dict = {}
    pop_dict['CHROM:POS'] = '{ALLELE:FREQ}'
    n_sites = 0
    pop_sites = []
    with open(pop_frq) as f:
        for line in f:
            n_sites +=1
            parsed_line = line.split('\t')
            pop_sites.append(':'.join([parsed_line[0], parsed_line[1]]))
            try:
                pop_dict[':'.join([parsed_line[0],
                                   parsed_line[1]
                                  ])] = [float(j[2:]) for j in [
                                                 parsed_line[4],
                                                 parsed_line[5].replace('\n','')
                                                               ]]
            except IndexError:
                pass
    
    # Initialize a table to store p-values.
    
    p_vals = []
    
    # For each of the pop_frq sites...
    
    for i in list(range(1,n_sites)):
         
        # Extract the column representing the current site from the
        # empirical distribution file and parse it into numerical data.
        
        dist = []
        len_dist = 0
        with open(dist_path) as csvfile:
            dist_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            next(dist_reader)
            for row in dist_reader:
                len_dist +=1
                try:
                    dist.append([float(j[2:]) for j in row[i].split(' ')])
                except ValueError:
                    dist.append(2*[float('NaN')])
        
        # Convert the extracted data to chi-squared statistics, using the 
        # grand population frequency data as an expectation.
        
        for j in list(range(len_dist)):
            observed = numpy.array(dist[j])
            expected = numpy.array(all_dict[pop_sites[i]])
            chisq = stats.chisquare(observed, expected)[0]
            dist[j] = chisq
        
        # Compute the same statistic for the frequency data computed
        # on the input population.
        
        pop_obs = numpy.array(pop_dict[pop_sites[i]])
        pop_exp = numpy.array(all_dict[pop_sites[i]])
        pop_chisq = stats.chisquare(pop_obs, pop_exp)[0]
        
        # Compute a p-value for the population frequency based on the empirical
        # distribution.

        p_val = sum([1 if val >= pop_chisq else 0 for val in dist])/len_dist
        
        # Append this result to the p_vals table
        
        p_vals.append([pop_sites[i], p_val])
    
    # Carry out a Benjamini-Hochberg procedure on the recorded p-values and
    # tag each site as 'POP-SPECIFIC' or 'INCONCLUSIVE'. Then write the results
    # to a tab-delimited output file.
    
    # Sort the sites by p-value
    
    p_vals.sort( key=lambda entry: entry[1])
    
    # Apply Benjamini-Hochberg
    
    cutoff = 0
    for i in list(range(n_sites-1)):
        if p_vals[i][1] <= i/n_sites * fdr:
            cutoff = int(i)
    for i in list(range(n_sites-1)):
        if i <= cutoff:
            p_vals[i].append(pop+'-SPECIFIC')
        else:
            p_vals[i].append('INCONCLUSIVE')
            
    # Sort the sites by coordinates.        
    
    p_vals.sort( key=lambda entry: entry[0])
    
    # Write the results to output.
    
    with open(verdict_path, 'w') as f:
        f.write('\t'.join(['SITE', 'P_VALUE','VERDICT']))
        for entry in p_vals:
            f.write('\t'.join([str(x) for x in entry])+'\n')
            
    # Print a summary of the results.
    
    print(''.join(['Analysis complete: discovered ', str(cutoff), ' ', pop, 
                   '-specific sites and ', str(n_sites-cutoff),
                   ' inconclusive sites.'
                  ]))
                  
                  
                  
                  
#------------------------------------------------------------------------------#
#                             10. Cleaning up                                  #
#------------------------------------------------------------------------------#       

def clean_up(inpath, pop):
    print('Cleaning up.')
    
    # Remove temporary files.
    
    os.system(' '.join(['rm', inpath[::-1][4:][::-1]+'.vcf',
                inpath[::-1][4:][::-1]+'_ALL.frq',
                inpath[::-1][4:][::-1]+'_ALL.log',
                inpath[::-1][4:][::-1]+'_ALL.vcf',
                inpath[::-1][4:][::-1]+'_ALL_sort.vcf',
                inpath[::-1][4:][::-1]+'_ALL_sort.vcf.vcfidx',
                inpath[::-1][4:][::-1]+'_dist.csv',
                inpath[::-1][4:][::-1]+'_'+pop+'.frq',
                inpath[::-1][4:][::-1]+'_'+pop+'.log',
                inpath[::-1][4:][::-1]+'_'+pop+'.vcf',
                inpath[::-1][4:][::-1]+'_'+pop+'_sort.vcf',
                inpath[::-1][4:][::-1]+'_'+pop+'_sort.vcf.vcfidx',
                inpath[::-1][4:][::-1]+'_rand.frq',
                inpath[::-1][4:][::-1]+'_rand.log',
                inpath[::-1][4:][::-1]+'_rand.vcf',
                inpath[::-1][4:][::-1]+'_rand.vcf.vcfidx',
                inpath[::-1][4:][::-1]+'_tabix.txt',
                '{}/{}{}{}'.format(os.getcwd(), '20100804.', 'ALL', '.panel'),
                '{}/{}{}{}'.format(os.getcwd(), '20100804.', pop, '.panel'),
                '/'.join([os.getcwd(),
                          'ALL.2of4intersection.20100804.genotypes.vcf.gz.tbi'
                         ])
                       ]))
       
       
#------------------------------------------------------------------------------#
#                                   Main                                       #
#------------------------------------------------------------------------------#

# Extract a specific population from the population panel.

subset_pop(args.pop)

# Decide between sites and genes and convert the input to tabix input format.

if args.gtf == None:
    convert_sites(args.inpath, args.buff)
else:
    convert_genes(args.inpath, args.gtf, args.buff)
    
# Execute tabix, generating a horizontally sliced .vcf file.

execute_tabix(args.inpath)

# Generate a population-subsetted (vertically sliced) .vcf file.

subset_vcf(args.inpath, args.pop)

# Compute allele frequencies on the full (horizontally sliced) and subsetted 
# (horizontally and vertically sliced) .vcf files.

compute_frequencies(args.inpath, args.pop)

# Identify population specific sites via chi-squared test or permutation test

if args.chisq:
    chi_squared_test(args.fdr, args.inpath, args.pop)
else:
    empirical_dist(args.inpath, args.pop, args.rep)
    permutation_test(args.fdr, args.inpath, args.pop, args.rep)
    
# Remove temporary files.

clean_up(args.inpath, args.pop)
