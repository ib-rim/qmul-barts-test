import os

#This Python file path
file_path = os.path.abspath("test.ipynb") 

#TSV file path
gc = os.path.join(os.path.dirname(file_path), "./gene_coordinates.tsv")
#VCF file path
gv = os.path.join(os.path.dirname(file_path), "./test.vcf")

#Add headings to output file
outputFile = open("output.tsv", "w")
outputFile.write('#CHROM\tPOS\tREF\tALT\tINFO\tGENE')

#Read genetic variations file
with open(gv) as gv_file:
    chroms = {}

    for gv_line in gv_file:
        genes = []
        genesString = ''
        #Skip headings/irrelevant rows
        if(not gv_line.startswith('#')):
            gv_values = gv_line.split('\t')
            gv_chrom = gv_values[0]
            gv_pos = int(gv_values[1])
            gv_ref = gv_values[3]
            gv_alt = gv_values[4]
            gv_info = gv_values[7]

            #Only loop through coordinates with matching chroms to check pos
            if gv_chrom in chroms:
                for gc_values in chroms[gv_chrom]:
                    gc_gene = gc_values[0]
                    gc_start_pos = int(gc_values[2])
                    gc_end_pos = int(gc_values[3])
                    
                    #If between start and end pos
                    if gv_pos >= gc_start_pos and gv_pos <= gc_end_pos:
                        #Add to gene list
                        genes.append(gc_gene)
            else:
                #Loop through all lines of genomic coordinates file once to save gene according to chrom
                with open(gc) as gc_file:
                    for gc_line in gc_file:

                        #Skip headings/irrelevant rows
                        if(not gc_line.startswith('#')):
                            gc_values = gc_line.split('\t')
                            gc_gene = gc_values[0]
                            gc_chrom = gc_values[1]
                            gc_start_pos = int(gc_values[2])
                            gc_end_pos = int(gc_values[3])

                            #Find matching CHROMs
                            if gc_chrom == gv_chrom:
                                #If between start and end pos
                                if gv_pos >= gc_start_pos and gv_pos <= gc_end_pos:
                                    #Add to gene list
                                    genes.append(gc_gene)

                            #Add to chroms dict      
                            if gc_chrom not in chroms:
                                chroms[gc_chrom] = []
                            chroms[gc_chrom].append(gc_values)

            #Convert list of genese to string
            #with format: 'gene1;gene2;gene3' or '.' if empty
            genesString = ';'.join(genes)
            if genesString == '':
                genesString = '.'
                
            #Add to output file in format:  #CHROM POS REF ALT INFO GENE
            outputFile.write(f'\n{gv_chrom}\t{gv_pos}\t{gv_ref}\t{gv_alt}\t{gv_info}\t{genesString}')

outputFile.close()
