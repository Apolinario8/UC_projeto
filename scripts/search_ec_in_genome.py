from genome_annotation import FunctionalAnnotation


def main():
    genome_annotation = FunctionalAnnotation(r"C:\Users\gonca\OneDrive\Documentos\GitHub\UC_projeto\DB")
    genome_annotation.identify_gene_by_homology_from_ec("blastp", "6.2.1.76")



if __name__ == '__main__':
    main()