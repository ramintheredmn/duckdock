import biotite.database.rcsb as rcsb
import redo
import pypdb

def get_pdb(uniprot_id, exp_method="X-RAY DIFFRACTION", max_res=3, min_ligand_w=100):

    query_by_uniprot = rcsb.FieldQuery("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", exact_match=uniprot_id)
    query_by_expmethpd = rcsb.FieldQuery("exptl.method", exact_match=exp_method)
    query_by_res = rcsb.FieldQuery("rcsb_entry_info.resolution_combined", less_or_equal=max_res)
    query_by_ligand_mw = rcsb.FieldQuery("chem_comp.formula_weight", molecular_definition=True, greater=min_ligand_w)

    query = rcsb.CompositeQuery(
        [
            query_by_uniprot,
            query_by_res,
            query_by_expmethpd,
            query_by_ligand_mw
        ],
        "and"
    )

    pdb_ids = rcsb.search(query)
    print(len(pdb_ids))
    print("PDB IDs")
    return pdb_ids

@redo.retriable(attempts=10, sleeptime=2)
def describe_one_pdb_id(pdb_id):
    """Fetch meta information from PDB."""
    described = pypdb.describe_pdb(pdb_id)
    if described is None:
        print(f"! Error while fetching {pdb_id}, retrying ...")
        raise ValueError(f"Could not fetch PDB id {pdb_id}")
    return described

uniprot_id = input("Enter the uniprot id of the protein of interest: ")
pdb_ids = get_pdb(uniprot_id)
pdb_descs = [describe_one_pdb_id(pdb_id) for pdb_id in pdb_ids]

for desc in pdb_descs[:3]:
    print(desc)
