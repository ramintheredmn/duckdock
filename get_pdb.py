import biotite.database.rcsb as rcsb
import redo
import pypdb
import pandas as pd
import inquirer

@redo.retriable(attempts=10, sleeptime=2)
def get_pdb(uniprot_id, exp_method="X-RAY DIFFRACTION", max_res=float(3), min_ligand_w=100, max_chain=None):

    query_by_uniprot = rcsb.FieldQuery("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", exact_match=uniprot_id)
    query_by_expmethpd = rcsb.FieldQuery("exptl.method", exact_match=exp_method)
    query_by_res = rcsb.FieldQuery("rcsb_entry_info.resolution_combined", less_or_equal=max_res)
    query_by_ligand_mw = rcsb.FieldQuery("chem_comp.formula_weight", molecular_definition=True, greater=min_ligand_w)
    query_by_numberOfChains = rcsb.FieldQuery("rcsb_entry_info.deposited_polymer_entity_instance_count", less_or_equal=max_chain)

    if max_chain:
        query = rcsb.CompositeQuery(
            [
                query_by_uniprot,
                query_by_res,
                query_by_expmethpd,
                query_by_ligand_mw,
                query_by_numberOfChains
            ],
            "and"
        )
    else:
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

def extract_pdb_info(pdb_metadata_list):
    results = []
    
    for metadata in pdb_metadata_list:
        pdb_id = metadata.get('rcsb_id', 'Unknown')  # Get PDB ID
        
        resolution = 'N/A'
        ligand_count = 0
        ligand_names = []
        
        if 'rcsb_entry_info' in metadata and 'resolution_combined' in metadata['rcsb_entry_info']:
            resolution = metadata['rcsb_entry_info']['resolution_combined'][0] if metadata['rcsb_entry_info']['resolution_combined'] else 'N/A'
        elif 'rcsb_entry_info' in metadata and 'diffrn_resolution_high' in metadata['rcsb_entry_info']:
            resolution = metadata['rcsb_entry_info']['diffrn_resolution_high'].get('value', 'N/A')
        
        if 'rcsb_entry_info' in metadata:
            ligand_count = metadata['rcsb_entry_info'].get('nonpolymer_entity_count', 0)
        
        if 'rcsb_binding_affinity' in metadata:
            ligand_names = list(set([entry['comp_id'] for entry in metadata['rcsb_binding_affinity']]))

        if 'rcsb_entry_info' in metadata:
            num_chains = metadata['rcsb_entry_info'].get('deposited_polymer_entity_instance_count', 0)
        else:
            num_chains = None

        results.append({
            'pdb_id': pdb_id,
            'resolution': resolution,
            'ligand_count': ligand_count,
            'ligand_names': ', '.join(ligand_names) if ligand_names else 'None',
            'num_chains': num_chains
        })
    
    return results

max_res = input("Enter the max resolution you want: (empty for 3) ")
max_chain = input("Enter the max number of chains: (empty for not limiting) ")
uniprot_id = input("Enter the uniprot id of the protein of interest: ")

if len(max_res)>=1:
    pdbs = get_pdb(uniprot_id, max_res=float(max_res))
elif len(max_res) >=1 and len(max_chain)>=1:
    pdbs = get_pdb(uniprot_id, max_res=float(max_res), max_chain=int(max_chain))
elif len(max_chain) >= 1:
    pdbs = get_pdb(uniprot_id, max_chain=float(max_chain))
else:
    pdbs = get_pdb(uniprot_id)

pdb_descs = [describe_one_pdb_id(pdb_id) for pdb_id in pdbs]
res = extract_pdb_info(pdb_descs)
df = pd.DataFrame.from_dict(res)
df_m = df.query("ligand_names != 'None' and ligand_count >= 1").reset_index(drop=True)

choises = [f"{i[0]}) {i[1]['pdb_id']}, resulotion: {i[1]['resolution']}, Chain count: {i[1]['num_chains']}, Ligand: {i[1]['ligand_names']}" for i in df_m.iterrows()]

questions = [inquirer.List(
    'row',
    message='Select a protein',
    choices=choises,
)]

answer= inquirer.prompt(questions)

print(answer.values()[0].split(')')[0])

