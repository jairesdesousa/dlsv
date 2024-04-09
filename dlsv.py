import rdkit
from rdkit import Chem
import copy
from cddd.inference import InferenceModel

inference_model = InferenceModel()

r = Chem.SmilesMolSupplier("smiles.smi", nameColumn=-1, titleLine=False)
f = open ("dlsv.csv","w")

for m in r:
    smiles_embedding_m = inference_model.seq_to_emb(Chem.MolToSmiles(m))
    for atomo in m.GetAtoms():
        if atomo.GetSymbol() == 'H' or atomo.GetSymbol() == 'P':
            continue
        mnew = copy.copy(m)
        idatom = atomo.GetIdx()
        atomtochange = mnew.GetAtomWithIdx(idatom)
        atomtochange.SetAtomicNum(15)
        mnew_smiles = Chem.MolToSmiles(mnew)
        smiles_embedding_mnew = inference_model.seq_to_emb(mnew_smiles)
        for i in smiles_embedding_mnew-smiles_embedding_m:
            for j in i:
                f.write(f'{j},')
        f.write (mnew_smiles + ',' 
                  + str(atomo.GetSymbol())+','
                  + str(atomo.GetSymbol())+ str(atomo.GetTotalDegree())+'H'+str(atomo.GetTotalNumHs())+','
                  + str(atomo.GetSymbol())+ str(atomo.GetHybridization())+','
                  + str(atomo.GetIsAromatic())+'\n')
f.close

