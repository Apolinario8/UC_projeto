import pandas as pd
###  show all columns of a dataframe
pd.set_option('display.max_columns', None)




def main():
    transit = parse_results(r"C:\Users\gonca\OneDrive\Ambiente de Trabalho\Merlin\merlin-4.0.12\ws\cvulgaris\3077\transyt\results\scoresMethod1.txt")
    res = count_family(transit)
    print(transit)
    print(res)

def load_transyt_results(filenames):
    results , method2 = {}, {}
    for filename in filenames:
        if filename.endswith('scoresMethod1.txt'):
            results = parse_results(filename)
        elif filename.endswith('scoresMethod2.txt'):
            method2 = parse_method2(filename)
        else:
            raise Exception('Invalid file name')
    results.update(method2)
    return results

def count_family(results):
    res = {}
    for key, value in results.items():
        family = value.split(".")[0]
        if family in res.keys():
            res[family] += 1
        else:
            res[family] = 1
    return res

def load_transportdb_results(filename):
    transportdb_results = pd.read_csv(r"C:\Users\Bisbii\Desktop\TranSyt\case_study\L0lFMTn9.csv", header=None, index_col=0)
    transportdb_results.columns = ["substrate", "col1", "col2", "col3", "col4", "family"]
    transportdb_results = transportdb_results[(~transportdb_results["family"].str.startswith("9"))]
    transportdb_results_temp = transportdb_results.drop(columns=["substrate", "col1", "col2", "col3", "col4"], axis=1)
    ## discard rows where family starts by 9
    transportdb_results_dict = transportdb_results_temp.to_dict(orient="index")
    for key, value in transportdb_results_dict.items():
        transportdb_results_dict[key] = value["family"]
    return transportdb_results, transportdb_results_dict


def parse_results(file_path):
    results = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            key = lines[i].strip()[1:]
            value = lines[i + 1].strip().split('-')[0]
            results[key] = value
            i += 2
        else:
            i += 1
    return results

def parse_method2(file_path):
    genes = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            genes[line.split("\t")[0]] = ""
    return genes

if __name__ == '__main__':
    main()