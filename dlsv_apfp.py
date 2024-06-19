import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
import copy
import pandas as pd
import numpy as np
from cddd.inference import InferenceModel
from cddd.preprocessing import preprocess_smiles

size=512
apgen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=size, maxDistance=2)

inference_model = InferenceModel()
r = Chem.SmilesMolSupplier("smiles.smi", nameColumn=-1, titleLine=False)
fdlsv = open ("dlsv.csv","w")
fapfp = open ("apfp.csv","w")

for m in r:
    AllChem.ComputeGasteigerCharges(m)
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
                fdlsv.write(f'{j},')
        gcharge=str(atomo.GetProp('_GasteigerCharge'))
        fdlsv.write (mnew_smiles + ',' 
                    + str(atomo.GetSymbol())+ str(atomo.GetTotalDegree())+'H'+str(atomo.GetTotalNumHs())+','
                    + str(atomo.GetIdx())+ ','
                    + gcharge +'\n')
        apfp=apgen.GetFingerprint(m,fromAtoms=[atomo.GetIdx()])
        apfp=np.array(apfp)
        for i in range (size):
            fapfp.write(f'{apfp[i]},')
        fapfp.write(f'{mnew_smiles},{atomo.GetSymbol()},{atomo.GetIdx()},{gcharge}\n')
fdlsv.close
fapfp.close


