from genome_annotation import FunctionalAnnotation


def main():
    genome_annotation = FunctionalAnnotation(r"C:\Users\gonca\Documents\GitHub\UC_projeto\DB") 
    genome_annotation.identify_gene_by_homology_from_ec("blastp", "2.4.2.1")



if __name__ == '__main__':
    main()