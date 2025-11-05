# PCB (Protein Charge Blockiness)
The library is designed to determining protein charge blockiness.

Liquid-liquid phase separation (LLPS) is a process in which biomolecules (primarily proteins and nucleic acids) in a cell spontaneously separate from the surrounding solution, forming biomolecular condensates that function as non-membrane organelles (e.g., the nucleolus, stress granules, and nuclear Cajal bodies).

Weak interactions are the key driving force of LLPS. Although hydrophobic and π-π interactions play a role, electrostatic interactions make the largest contribution.

Charge distribution patterns in the linear amino acid sequence of a protein determine its ability to undergo phase separation and its choice of interaction partners.

## Installing
1. Cloning repository
```bash 
git clone https://github.com/BioInfDev/PCB.git
```
2. Change current directory to pcb
```bash
cd ./PCB
```
3. Create and activate venv
```bash
sudo apt install python3.12-venv
python3 -m venv venv
source venv/bin/activate
```  
4. Installing dependencies
```bash
pip install -r requirements.txt
```
5. Run the example
```bash
./example.py
```

## Example
The example analyzes the charge blocking of the nucleolar protein Treacle (Uniprot ID: Q14328).
<img width="1920" height="1080" alt="Screenshot from 2025-11-05 15-04-17" src="https://github.com/user-attachments/assets/e13d2cc8-7f80-447a-95ac-aab63fee8124" />
