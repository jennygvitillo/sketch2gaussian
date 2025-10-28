# sketch2gaussian - molecular structure converter - From ChemDraw to Gaussian input files

**Python** tool for converting 2D molecular structures from ChemDraw to professional Gaussian input files with accurate 3D geometries. From simple drawings to quantum-ready calculations.

## Features

- **2D to 3D Conversion**: Transform ChemDraw (.mol) files into Gaussian (.gjf) input files
- **Automatic Hydrogen Addition**: Adds missing hydrogen atoms with correct valences  
- **Accurate 3D Geometries**: Uses RDKit's ETKDG algorithm for chemically correct bond lengths and angles
- **SMILES Integration**: Includes SMILES strings in Gaussian files for molecular traceability
- **Batch Processing**: Processes multiple .mol files in a folder automatically
- **Interactive GUI**: User-friendly folder selection dialog

## Input Requirements
- ChemDraw molecular structure files (.mol format)
- Open Babel installation (for SMILES conversion)
- RDKit Python package (for 3D coordinate generation)

## How to run it
(on Windows, for dummies) Detailed instructions on how to run the script are reported at the beginning of the script.
for the other users: like all the other Python tools.

## Output files
**Gaussian Input Files (.gjf)**

**Features included in output files:**
- **Professional 3D Coordinates**: Chemically accurate geometries with proper bond lengths
- **SMILES Integration**: Each file includes the SMILES string for molecular reference
- **Standard Computational Parameters**: Pre-configured with B3LYP/def2-TZVP and common options
- **Automatic Hydrogen Addition**: Complete molecular structures ready for calculation

 ``` input_folder/
├── molecule1.mol
├── molecule2.mol
└── gjf_files_rdkit/ # Created automatically
├── molecule1.gjf # Gaussian input with 3D coordinates
└── molecule2.gjf # Includes SMILES: C=Cc1ccccc1
 ``` 

## Citation
J.G. Vitillo, 'sketch2gaussian', https://github.com/jennygvitillo/sketch2gaussian (2025).

## Italiano
Strumento **Python** per convertire strutture molecolari 2D da ChemDraw in file di input Gaussian professionali con geometrie 3D accurate. Da disegni semplici a calcoli quantistici pronti all'uso.

**Requisiti di Input**
- File di strutture molecolari ChemDraw (formato .mol)
- Installazione di Open Babel (per conversione SMILES)
- Pacchetto Python RDKit (per generazione coordinate 3D)

- **Conversione 2D a 3D** - trasforma file ChemDraw (.mol) in file di input Gaussian (.gjf)
- **Aggiunta Automatica Idrogeni** - aggiunge atomi di idrogeni mancanti con valenze corrette
- **Geometrie 3D Accurate** - utilizza l'algoritmo ETKDG di RDKit per lunghezze di legame e angoli chimicamente corretti
- **Integrazione SMILES** - include stringhe SMILES nei file Gaussian per tracciabilità molecolare
- **Elaborazione in Batch** - processa automaticamente multipli file .mol in una cartella
- **Interfaccia Grafica Interattiva** - dialogo user-friendly per selezione cartelle

**Dipendenze**: rdkit, openbabel
**Istruzioni dettagliate** su come usare lo script sono riportate nella parte iniziale dello script stesso (file: sketch2gaussian.py).
