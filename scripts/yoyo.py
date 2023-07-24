from cobra import Model
from cobra.manipulation import remove_genes


def get_reactions(id_lipid, reactions_list, visited):
    lipid = model_cyano.metabolites.get_by_id(id_lipid + "__cytop")
    for reaction in lipid.reactions:
        if lipid in reaction.products or reaction.lower_bound < 0:
            reactions_list.add(reaction)
            for reactant in reaction.reactants:
                if reactant.id not in visited and (reactant.id.startswith("BMGC") or reactant.id.startswith("mgdg") or reactant.id.startswith("dgdg") or reactant.id.startswith("sqdg")):
                    visited.add(reactant.id)
                    reactions_to_add = get_reactions(reactant.id.split("__")[0], reactions_list, visited)
                    for reaction_to_add in reactions_to_add:
                        reactions_list.add(reaction_to_add)
    return reactions_list


from cobra.io import read_sbml_model, write_sbml_model

model_cyano = read_sbml_model('model_lipids_cyano.xml')

id_lipid = ["BMGC31499", "pg18111Z1613E", "BMGC456", "BMGC765976", "BMGC7150", "BMGC765977", "BMGC3716", "BMGC765978", #PG
            
            "pe1801835Z9Z12Z", "pe1801845Z9Z12Z15Z", "pe18111Z1835Z9Z12Z", "pe18111Z1845Z9Z12Z15Z", "pe1819Z1835Z9Z12Z", "pe1819Z1845Z9Z12Z15Z", "pe1829Z12Z1835Z9Z12Z", #PE
            "BMGC4854", "BMGC42965", #PI

            "pc1829z12z1835z9z12z", "BMGC415", "pc1819z1819z", "BMGC302", "BMGC45682", "pc1801835z9z12z", "BMGC10694", "BMGC10692", "pc18111z1835z9z12z", #PC

            "tag16018111Z160", "BMGC98268", "tag16018111Z18111Z", "BMGC80148", "tag16018111Z1835Z9Z12Z", "tag16018111Z1845Z9Z12Z15Z", "tag1601819Z160", "BMGC98270", 
            "BMGC80154", "tag1601819Z1819Z", "tag1601819Z1835Z9Z12Z", "tag1601819Z1845Z9Z12Z15Z", "BMGC102109", "tag1801819Z180", "BMGC111412", "tag1801819Z1819Z", 
            "tag1801819Z1835Z9Z12Z", "tag1801819Z1845Z9Z12Z15Z", "BMGC208464", "BMGC208461", "tag18111Z18111Z18111Z", "tag18111Z18111Z1819Z", "tag18111Z18111Z1835Z9Z12Z",
            "tag18111Z18111Z1845Z9Z12Z15Z", "tag18111Z1819Z160", "BMGC208463", "tag18111Z1819Z18111Z", "tag18111Z1819Z1819Z", "tag18111Z1819Z1835Z9Z12Z", 
            "tag18111Z1819Z1845Z9Z12Z15Z", "tag1819Z18111Z160", "BMGC111562", "tag1819Z18111Z18111Z", "BMGC76658", "tag1819Z18111Z1835Z9Z12Z", "tag1819Z18111Z1845Z9Z12Z15Z",
            "tag1819Z1819Z160", "tag1819Z1819Z180", "tag1819Z1819Z18111Z", "tag1819Z1819Z1819Z", "tag1819Z1819Z1835Z9Z12Z", "tag1819Z1819Z1845Z9Z12Z15Z", #TAG

            "sqdg160", "sqdg18111Z160", "sqdg1819Z160", "sqdg1829Z12Z160", "sqdg1839Z12Z15Z160", #SQDG

            "mgdg1819Z160", "mgdg1819Z1617Z", "mgdg1819Z1619Z", "mgdg1819Z1627Z10Z", "mgdg1819Z1637Z10Z13Z", "mgdg1829Z12Z160", "mgdg1829Z12Z1617Z", "mgdg1829Z12Z1619Z",
            "mgdg1829Z12Z1627Z10Z", "mgdg1829Z12Z1634Z7Z10Z", "mgdg1829Z12Z1637Z10Z13Z", "mgdg1829Z12Z1644Z7Z10Z13Z", "mgdg1839Z12Z15Z160", "mgdg1839Z12Z15Z1627Z10Z",
            "mgdg1839Z12Z15Z1634Z7Z10Z", "mgdg1839Z12Z15Z1637Z10Z13Z", "mgdg1839Z12Z15Z1644Z7Z10Z13Z", #MGDG

            "dgdg1819Z160", "dgdg1819Z1617Z", "dgdg1819Z1619Z", "dgdg1819Z1627Z10Z", "dgdg1819Z1634Z7Z10Z", "dgdg1819Z1637Z10Z13Z", "dgdg1829Z12Z160", "dgdg1829Z12Z1617Z",
            "dgdg1829Z12Z1619Z", "dgdg1829Z12Z1627Z10Z", "dgdg1829Z12Z1634Z7Z10Z", "dgdg1829Z12Z1637Z10Z13Z", "dgdg1839Z12Z15Z160", "dgdg1839Z12Z15Z1627Z10Z", 
            "dgdg1839Z12Z15Z1634Z7Z10Z", "dgdg1839Z12Z15Z1637Z10Z13Z" #DGDG
            ]

id_metabolites = [metabolite.id for metabolite in model_cyano.metabolites]
results = set()
not_in_model = set()

for lipid in id_lipid:
    if lipid + "__cytop" in id_metabolites:
        reactions = get_reactions(lipid, set(), set())
        for reaction in reactions:
            for met in reaction.metabolites:
                met.compartment = "cytop"
            results.add(reaction)
    else:
        print(lipid + "__cytop" + " not in model")
        not_in_model.add(lipid)

with open("not_in_model.txt", "w") as f:
    for lipid in not_in_model:
        f.write(lipid + "\n")

new_model = Model()

new_model.add_reactions(results)

groups_to_add = []
for group in model_cyano.groups:
    if any([reaction in results for reaction in group.members]):
        to_remove = []
        for reaction in group.members:
            if reaction not in results:
                to_remove.append(reaction)
        group.remove_members(to_remove)
        groups_to_add.append(group)

new_model.add_groups(groups_to_add)
# write

remove_genes(new_model, new_model.genes, remove_reactions=False)

write_sbml_model(new_model, "model_lipids_cyano_yo.xml")
