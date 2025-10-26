# sketch2gaussian
This script converts 2D molecular structures from ChemDraw (.mol files) into 3D structures with realistic bond lengths and geometries. I have chosen as output format input files for Gaussian 16 (.gjf) because it is a largely used program for quantum mechanics calculations but the XYZ coordinates can be very easily used for generate manually a XYZ file by cut and paste. It adds all the hydrogens!! It uses a hybrid approach combining Open Babel for SMILES conversion and RDKit's ETKDG algorithm for accurate 3D coordinate generation. The instructions for Windows users and that have not used Python before are reported at the beginning of the file. Final structure of the folder:

your_folder/

├── molecule1.mol

├── molecule2.mol

└── gjf_files_rdkit/ # Created automatically

     ├── molecule1.gjf         # Gaussian input

     └── molecule2.gjf
This project is licensed under the MIT License (see the LICENSE file for details).

Questo script converte strutture molecolari 2D provenienti da ChemDraw (file .mol) in strutture 3D con lunghezze di legame e geometrie realistiche. Ho scelto come formato di output i file di input per Gaussian 16 (.gjf), poiché è un programma ampiamente utilizzato per i calcoli di meccanica quantistica. Tuttavia, le coordinate XYZ generate possono essere facilmente riutilizzate per creare manualmente un file .xyz tramite copia e incolla. Aggiunge automaticamente tutti gli atomi di idrogeno, quindi le strutture possono avere gli idrogeni impliciti! Utilizza un approccio ibrido che combina Open Babel per la conversione in SMILES e l’algoritmo ETKDG di RDKit per la generazione accurata delle coordinate 3D. Le istruzioni per utenti Windows alle prime armi con Python sono riportate all’inizio del file. Struttura finale della cartella:

cartella_originaria/

├── molecule1.mol

├── molecule2.mol

└── gjf_files_rdkit/ # Created automatically

  ├── molecule1.gjf         # Gaussian input

  └── molecule2.gjf
Questo progetto ha licenza MIT (vedi il file LICENSE per dettagli).

**Cite it as:
J. G. Vitillo, "sketch2gaussian", https://github.com/jennygvitillo/sketch2gaussian (2025).**
